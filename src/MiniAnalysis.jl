module MiniAnalysis
using Base
using LsqFit
using Statistics
using Printf
using FindPeaks1D#
using Random, Distributions
using DSP
using ArraysOfArrays
using ElasticArrays

using Distributed
using SharedArrays  # important for parallel/looping
using ProgressMeter
using Base.Threads
using FFTW
using ForwardDiff
using OnlineStats
#using BenchmarkTools  Throws error in package

#ENV["MPLBACKEND"] = "MacOSX"
using Plots
pyplot()


export find_events
export testdetector, testtemplate, make_test_waveforms

function clements_bekkers(data::Array{Float64, 1}, pars::NamedTuple{(:t, :dt), Tuple{Vector{Float64}, Float64}}) # template::Array{Float64, 1}, tau, dt)
    # , data2::Array{Float64, 1}, u::Array{Int64, 1}, sse, scale, offset, detcrit)
    #=
    cb algorithm
    =#
    ## Prepare arrays
    # println(size(data), size(template))
    # extract from structure
    template = pars.t
    # other values are not needed
    n_template = size(template)[1]
    n_data = size(data)[1]
    n_dt = n_data - n_template
    sume::Float64 =  sum(template)
    sume2::Float64 = sum(template .^ 2)
    # data = data.*1e12
    data_2 = Array{Float64, 1}(undef, size(data)[1])
    data_2 = data .^ 2
    sumy::Float64 = sum(data[1:n_template])
    sumey::Float64 = 0.
    s::Float64 = 0.
    c::Float64 = 0.
    u = Array{Int64, 1}(undef, (n_template))
    j::Int64=0
    sse = Array{Float64, 1}(undef, (n_template))
    sumy2::Float64 = sum(data_2[1:n_template])
    scale = zeros(Float64, (n_data, 1))
    offset = zeros(Float64, (n_data, 1))
    detcrit = zeros(Float64, (n_data, 1))
    for i = 1:n_dt
        k = i
        if i > 1
            sumy = sumy .+ data[k+n_template] .- data[k]
            sumy2 = sumy2 .+ data_2[k+n_template] .- data_2[k]
        end
        m = k + n_template - 1
        sumey = sum(data[i:m] .* template)
        s =
            (sumey .- (sume .* sumy) ./ n_template) /
            (sume2 .- (sume .* sume) ./ n_template)
        c = (sumy .- (s .* sume)) ./ n_template
        # fitted_template = template * s + c
        sse = sumy2 + (s .* s .* sume2) + (n_template .* c .* c)
        sse -= 2.0 .* ((s .* sumey) + (c .* sumy) .- (s .* c .* sume))
        u = findall(sse .<= 0)
        if u != Int64[]  # prevent negative value...
            j = u[1]
            sse = -sse
        end
        # sum((data[i : i + n_template-1] .- fitted_template) .^ 2)
        detcrit[i] = s ./ sqrt(sse ./ (n_template - 1))
        scale[i] = s
        offset[i] = c
    end
    return scale, detcrit
end


function aj_deconvolve(data::Array{Float64, 1}, 
        pars::NamedTuple{(:template, :tau, :dt), Tuple{Vector{Float64}, Vector{Float64}, Float64}})
    #=
    pars is a tuple with template::Array{Float64, 1}, tau, dt
    =#
    template = pars.template
    tau = pars.tau
    dt = pars.dt
    
    llambda::Float64 = 5.0
    data = data .- mean(data)
    
    # Weiner filtering/convolution
    if size(template)[1] < size(data)[1]
        template = vcat(template, zeros(size(data)[1] - size(template)[1]))
    end
    H = fft(template)
    # if size(H)[1] < size(data)[1]
    #     H = vcat(H, zeros(size(data)[1] - size(H)[1]))
    # end

    outsize = size(data)[1]# - size(template)[1]
    quot = ifft(
        fft(data) .* conj(H) ./ (H .* conj(H) .+ llambda ^ 2.0)
    )
    detcrit = real(quot).*llambda
    return real(quot[1:outsize]), detcrit[1:outsize]
end


function rs_deconvolve(data::Array{Float64, 1}, pars::NamedTuple{(:tau, :dt), Tuple{Float64, Float64}})
    #pars is a tuple with: template, tau, dt)
    dt = pars.dt
    tau = pars.tau
    dVdt = diff(data)/dt
    dVdt = vcat(dVdt, dVdt[end])  # add a last point same as the first
    d = tau .* dVdt .+ data
    return dVdt, -d 

end

function measure_noise(data::Array{Float64, 1};  threshold=2.0, iterations=2)
    ## Determine the base level of noise

    if iterations > 1
        med = median(data)
        stdev = std(data)
        thresh = stdev * threshold
        datam = data[(data .> med) .| (data .< -med)]
        return measureNoise(datam, threshold, iterations-1)
    else
        return std(data)
    end
end
        
function zero_crossings(data::Array{Float64, 1}, 
    pars::NamedTuple{(:sign, :minLength, :minPeak, :minSum, :noiseThreshold),
            Tuple{Int64, Int64, Float64, Float64, Float64}})
    #=
    locate events of any shape in a signal. Works by finding regions of the signal
    that deviate from noise, using the area, duration and minimum peak beneath the deviation as the detection criteria.

    Makes the following assumptions about the signal:
      - noise is gaussian
      - baseline is centered at 0 (high-pass filtering may be required to achieve this).
      - no 0 crossings within an event due to noise (low-pass filtering may be required to achieve this)
      - Events last more than minLength samples
      Return an array of events where each row is (start, length, sum, peak)
   =#
    minLength = pars.minLength
    minPeak = pars.minPeak
    minSum = pars.minSum
    noiseThreshold = pars.noiseThreshold
    data = data .* convert(Float64, pars.sign)  # flip sign if needed
    ncrosses = findall((data[1:end-1] .> 0) .& (data[2:end] .<= 0))  ## get indices of points crossing 0 in the right direction
    ncrosses .+= 1 # adjust array indexing
    pcrosses = findall((data[1:end-1] .<= 0) .& (data[2:end] .> 0))  ## get indices of points crossing 0 in the right direction

    x = collect(range(0.0, size(data)[1]-1, step=1.0))

    # handle pairwise conditions for positive and negative crossings
    if ncrosses[1] < pcrosses[1]
        np = true
    else
        np = false
    end
    
    ev_index = zeros(Int64, size(ncrosses)[1])
    ev_amp = zeros(Float64, size(ncrosses)[1])
    ev_peak_index = zeros(Int64, size(ncrosses)[1])
    ev_sum = zeros(Float64, size(ncrosses)[1])
    k = 0
    l_n = size(ncrosses)[1]
    p_n = size(pcrosses)[1]
    println(l_n, " ", p_n)
    # if noisethreshold is > 0, then
    noise_value = 0
    # Fit gaussian to peak in size histogram, use fit sigma as criteria for noise rejection
    # stdev = measureNoise(data1)
    # #p.mark('measureNoise')
    # hist = np.histogram(events['sum'], bins=100)
    # #p.mark('histogram')
    # histx = 0.5*(hist[1][1:] + hist[1][:-1]) ## get x values from middle of histogram bins
    # #p.mark('histx')
    # fit = fitGaussian(histx, hist[0], [hist[0].max(), 0, stdev*3.0, 0])
    # #p.mark('fit')
    # sigma = fit[0][2]
    # if np.isnan(sigma):
    #     sigma = np.std(hist[0])
    # minSize = sigma * noiseThreshold
    
    #
    # Filter events on size, length and/or charge
    #
    for i = 1:l_n
        if np
            dwin = data[ncrosses[i]:pcrosses[i]]  # aligned
        else
            if i > (p_n-1)
                # println("I exceeds p_n")
                continue  # skip as no end in sight
            end
            dwin = data[ncrosses[i]:pcrosses[i+1]]  # always go from n to p
        end
        if size(dwin)[1] < minLength
            # println("rejected on length")
            continue
        end
        # println(minimum(dwin), " ", pars.minPeak)
        if (pars.minPeak != 0.0) & (maximum(dwin) < pars.minPeak)
            # println("rejected on min peak")
            continue
        end
        if (minSum > 0.) & (sum(dwin) < minSum)
            # println("rejected on sum ", minSum, " ", sum(dwin))
            continue
        end
        if (noiseThreshold > 0.) & (maximum(dwin) < noise_value)
            continue
        end
        # should deal with noise threshold here
        k += 1
        ev_index[k] = ncrosses[i]
        ev_amp[k], pki = findmin(dwin)
        ev_peak_index[k] = convert(Int64, pki)+ ncrosses[i]
        ev_sum[k] = sum(dwin)
        ev_peak_index[i] += ev_index[k] - 1
    end
    ev_index = ev_index[1:k]
    ev_amp = ev_amp[1:k]
    ev_sum = ev_sum[1:k]
    ev_peak_index = ev_peak_index[1:k]

    p1 = plot(x, data, linewidth=0.5, color="black")
    p1 = plot!(p1, x[ncrosses], data[ncrosses], seriestype = :scatter, markercolor="blue", markerstrokewidth=0, markersize=4)
    p1 = plot!(p1, x[pcrosses], data[ncrosses], seriestype = :scatter, markercolor="red", markerstrokewidth=0, markersize=4)
    p1 = plot!(p1, x[ev_peak_index], data[ev_peak_index], seriestype = :scatter, markercolor="cyan", markerstrokewidth=0, markersize=4)

    title = plot(
        title = @sprintf("%s Test", mode),
        grid = false,
        showaxis = false,
        yticks = false,
        xticks = false,
        bottom_margin = -50Plots.px,
    )
    plot(title, p1,layout=grid(2, 1, heights = [0.1, 0.45]))  # note use of show=true here - critical!
    plot!(size = (600, 600), show=true)
    show()
    crit = zeros(Float64, size(data)[1])
    crit[ev_index] .= 1.0
    return data, crit

end
    
function make_template(;
    taus = [1e-3, 3e-3],
    sign = 1,
    dt_seconds = 50e-6,
    tmax = 0.010,
    risepower = 1,
    tdelay = 0.0,
)

    idelay = floor(Int, tdelay / dt_seconds)
    t_psc = collect(range(0.0, tmax + tdelay, step = dt_seconds)) # collect(0:dt_seconds:(tmax+tdelay))

    Aprime = (taus[2] / taus[1])^(taus[1] / (taus[1] - taus[2]))
    template = Array{Float64,1}(undef, size(t_psc)[1])
    tm = (
        (1.0 ./ Aprime) .* ((1.0 .- (exp.(-t_psc ./ taus[1]))) .^ risepower) .*
        exp.(-t_psc ./ taus[2])
    )

    template = tm / maximum(tm)
    # end
    if sign > 0
        template_amax = maximum(template)
    else
        template = -template
        template_amax = minimum(template)
    end
    return t_psc, template
end

function test_template()
    t_template, template = make_template()
    println(size(t_template), size(template))
    plot(t_template, template)
    plot!(size = (500, 500))
end

function filter(tdat, idat; lpf = 3000, hpf = 0)
    #=
    Digitally filter the data in dat with either a
    lowpass filter (lpf) or an hpf filter (hpf)
    Values are in Hz if tdat is in seconds.
    Returns filtered idat
    =#
    fs = 1 / mean(diff(tdat[:, 1]))  # get fs in Hz
    if hpf == 0
        responsetype = Lowpass(lpf; fs = fs)
    else
        responsetype = Bandpass(hpf, lpf; fs = fs)
    end
    designmethod = Butterworth(8)
    idat = filt(digitalfilter(responsetype, designmethod), idat)
    return idat
end

function make_test_waveforms(;N=1, rate=5.0, maxt=5.0, lpf=3000.0, taus=[1e-3, 3e-3], sign=-1)
    #=
    Make test waveforms, consisting of events with exp distribution,
    variable amplitudes, convolved with a template
    
    N is the number of waveforms to generate
    =#
    
    dt_seconds = 50e-6

    tmax = 5 * taus[2]
    
    t_template, template =
        make_template(taus = taus, tmax = tmax, dt_seconds = dt_seconds, sign = sign)
    t_wave = collect(range(0.0, maxt, step = dt_seconds))
    n_t_wave = size(t_wave)[1]
    # plot(t_template, template)
    # show()
    nevents = floor(Int, maxt * rate)  # expected events in time window
    t_psc = zeros(Float64, n_t_wave, N)
    w_psc = zeros(Float64, n_t_wave, N)
    
    for i = 1:N
        Random.seed!(42+i)
        d = Exponential(1.0 / rate)
        i_noise = Normal(1.0, 0.3) * 20e-12
        event_variance = rand(i_noise, nevents)
        evtx = cumsum(rand(d, nevents))
        noise = Normal(0, 1) * 0.25e-12
        noise_add = rand(noise, size(t_wave)[1])
        evta = zeros(Float64, (size(t_wave)[1], 1)) .+ noise_add
        k = 1
        for ev in evtx
            ix = floor(Int, ev / dt_seconds)
            if ix == 0
                continue
            end
            if ix < size(evta)[1]
                evta[ix, 1] = 1.0 * event_variance[k]
                k += 1
            end
        end
        ievt = findall(evta .> 0.0)
        wv = conv(evta, template)[1:n_t_wave]  # clip to standard length
        tmax = size(wv)[1] * dt_seconds
        t_psc[:, i] = collect(range(0.0, tmax, length = size(wv)[1]))
        w_psc[:, i] = filter(t_psc, wv, lpf = lpf)
    end
    return t_psc, w_psc, dt_seconds, size(w_psc)[1], template

end


function identify_events(wave, crit, thr, dt; thresh = 3.0, ev_win = 5e-3)
    #=
    find events crossing threshold - these are onsets
    Operates on a single trace
    =#    
    ev = findall((crit[1:end-1] .> thr)) # .& (diff(crit) .> 0))
    if ev == []  # no events found!
        return [], [], thr
    end
    ev = insert!(ev, 1, 0)
    evn = findall(diff(ev) .> 1).+1
    ev = ev[evn]
    swin = floor(Int, ev_win / dt)

    # find maximum of event
    pks = zeros(Int32, (size(ev)[1]))
    maxn = size(ev)[1]
    k = 1
    for iev in ev
        if iev + swin < size(wave)[1]
            ipk = argmin(wave[iev:(iev+swin)])
            pks[k] = ipk + iev
            k += 1
        else
            maxn -= 1  # event goes beyond end, so delete it
        end
    end
    return ev[1:maxn], pks[1:maxn], thr
end

    
function testdetector(mode="CB";parallel=false, N=1, thresh=3.0, tau1=1, tau2=3, sign=1)
    N_tests = N
    t_psc, idat, dt_seconds, nwave, template = make_test_waveforms(N=N_tests, taus=[tau1*1e-3, tau2*1e-3], sign=sign)
    n_traces = size(idat)[2]
    ntemplate = nwave - size(template)[1]
    if mode == "CB"
        method = clements_bekkers
        pars = (t=template, dt=dt_seconds)
        
    elseif mode == "AJ"
        method = aj_deconvolve
        pars = (template=template, tau=[tau1, tau2], dt=dt_seconds)
    elseif mode == "RS"
        method = rs_deconvolve
        pars = (tau=tau1, dt=dt_seconds)  # tau instead of template
    elseif mode == "ZC"
        method = zero_crossings
        pars = (sign=sign, minLength=floor(Int, 1e-3/dt_seconds), minPeak=5e-12, minSum=0.0, noiseThreshold=2.5)
        
    else
        println("Method not recognized: ", mode)
    end
    if parallel
        s = SharedArray{Float64,2}((nwave, n_traces))
        c = SharedArray{Float64,2}((nwave, n_traces))
        @time @threads for i = 1:n_traces
                s[:,i], c[:, i] = method(idat[:,i], pars)
        end
    
    else
        s = Array{Float64,2}(undef, (nwave, n_traces))
        c = Array{Float64,2}(undef, (nwave, n_traces))
        @time for i = 1:n_traces
            s[:,i], c[:, i] = method(idat[:,i], pars)
        end
    end
    println("Detection complete")# draw the threshold line
    thr = std(c) * thresh  # all traces go into the threshold
    println("thr: ", thr, " ", thresh)

    ev = Array{Array{Int64}}(undef, (n_traces))
    pks = Array{Array{Int64}}(undef, (n_traces))
    thr_value = 0
    for i = 1:n_traces
        eva, pksa, thr_value = identify_events(idat[:, i], c[:, i], thr,  dt_seconds)
        ev[i] = eva
        pks[i] = pksa
    end

    thrline = [thr_value, thr_value]
    thrx = [0, maximum(t_psc[1:nwave, 1])]
    p1  = ""
    p2 = ""
    p3 = ""
    p4 = ""
    for i = 1:N_tests
        if i == 1
            p1 = plot(
                t_psc[:, i],
                idat[:, i],
                #color = "black",
                palette = :Dark2_5,
                linewidth = 0.5,
                title = "EPSC",
                titlefontsize = 10,
                legend = false,
            )
        else
            p1 = plot!(
                p1, t_psc[:, i],
                idat[:, i],
                #color = "black",
                palette = :Dark2_5,
                linewidth = 0.5,
                # title = "EPSC",
                # titlefontsize = 10,
                # legend = false,
            )
        end
        p1 = plot!(
            p1, t_psc[ev[i], i],
            idat[ev[i], i],
            seriestype = :scatter,
            markercolor = "yellow",
            markersize = 2.5,
        )

        p1 = plot!(
            p1, t_psc[pks[i], i],
            idat[pks[i], i],
            seriestype = :scatter,
            markercolor = "red",
            markersize = 2.5,
        )
        if i == 1
            p2 = plot(
            t_psc[1:size(s)[1], i],
            s[:, i],
            linecolor = "cyan",
            title = "Scale",
            titlefontsize = 10,
            legend = false,
        )
        else
            p2 = plot!(
            p2, t_psc[1:size(s)[1], i],
            s[:, i],
            linecolor = "cyan",
            title = "Scale",
            titlefontsize = 10,
            legend = false,
            )
        end
        if i == 1
            p3 = plot(
                t_psc[1:size(c)[1], i],
                c[:, i],
                #linecolor = "red",
                palette = :Dark2_5,
                title = "Criteria",
                titlefontsize = 10,
                legend = false,
                )
                p3 = plot!(p3, thrx, thrline, linecolor = "green", legend = false)
        else
            p3 = plot!(
                p3, t_psc[1:size(c)[1], i],
                c[:, i],
                #linecolor = "red",
                palette = :Dark2_5,
                title = "Criteria",
                titlefontsize = 10,
                legend = false,
                )
            p3 = plot!(p3, thrx, thrline, linecolor = "green", legend = false)
        end
        if i == 1
            p4 = plot(
                t_psc[1:size(template)[1]],
                template,
                linecolor = "blue",
                title = "Template",
                titlefontsize = 10,
                legend = false,
            )
        else
            p4 = plot!(p4,
                t_psc[1:size(template)[1]],
                template,
                linecolor = "blue",
                title = "Template",
                titlefontsize = 10,
                legend = false,
            )
        
        end
    end
    title = plot(
        title = @sprintf("%s Test", mode),
        grid = false,
        showaxis = false,
        yticks = false,
        xticks = false,
        bottom_margin = -50Plots.px,
    )
    l = @layout([a{0.1h}; b; c; d; e])
    plot(
        title,
        p1,
        p2,
        p3,
        p4,
        # layout = l,
        layout = grid(5, 1, heights = [0.1, 0.25, 0.25, 0.25, 0.15]),
    ) #, 0.30, 0.30]))
    plot!(size = (600, 600))
    # Plots.savefig("CB_test.pdf")
end


end
