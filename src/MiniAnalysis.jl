module MiniAnalysis
using Base
using LsqFit
using Statistics
using Printf
# using FindPeaks1D#
using Random, Distributions
using DSP
using ArraysOfArrays
using ElasticArrays
using Crayons.Box
using Distributed
using SharedArrays  # important for parallel/looping
# using ProgressMeter
using Base.Threads
using FFTW
using ForwardDiff
using DataFrames
# using OnlineStats
#using BenchmarkTools  Throws error in package

#ENV["MPLBACKEND"] = "MacOSX"
using Plots
pyplot()

include("LSPSPlotting.jl")
include("Acq4Reader.jl")

# externally visible routines:
export testdetector, testtemplate, make_test_waveforms
export detect_events, label_events
export unpack_events, events_to_dataframe, extract_events

function clements_bekkers(
    data::Array{Float64,1},
    pars::NamedTuple{(:t, :dt),Tuple{Vector{Float64},Float64}},
) # template::Array{Float64, 1}, tau, dt)
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
    sume::Float64 = sum(template)
    sume2::Float64 = sum(template .^ 2)
    # data = data.*1e12
    data_2 = Array{Float64,1}(undef, size(data)[1])
    data_2 = data .^ 2
    sumy::Float64 = sum(data[1:n_template])
    sumey::Float64 = 0.0
    s::Float64 = 0.0
    c::Float64 = 0.0
    u = Array{Int64,1}(undef, (n_template))
    j::Int64 = 0
    sse = Array{Float64,1}(undef, (n_template))
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


function aj_deconvolve(
    data::Array{Float64,1},
    pars::NamedTuple{(:template, :tau, :dt),Tuple{Vector{Float64},Vector{Float64},Float64}},
)
    #=
    pars is a tuple with template::Array{Float64, 1}, tau, dt
    =#
    template = pars.template
    tau = pars.tau
    dt = pars.dt

    llambda::Float64 = 5.0
    data = data .- mean(data)
    # println("size template: ", size(template))
    # Weiner filtering/convolution
    if size(template)[1] < size(data)[1]
        template = vcat(template, zeros(size(data)[1] - size(template)[1]))
    end
    # println(size(data), " ", size(template))
    H = fft(template)
    # if size(H)[1] < size(data)[1]
    #     H = vcat(H, zeros(size(data)[1] - size(H)[1]))
    # end

    outsize = size(data)[1]# - size(template)[1]
    quot = ifft(fft(data) .* conj(H) ./ (H .* conj(H) .+ llambda^2.0))
    detcrit = real(quot) .* llambda
    return real(quot[1:outsize]), detcrit[1:outsize]
end


function rs_deconvolve(
    data::Array{Float64,1},
    pars::NamedTuple{(:tau, :dt),Tuple{Float64,Float64}},
)
    #pars is a tuple with: template, tau, dt)
    dt = pars.dt
    tau = pars.tau
    dVdt = diff(data) / dt
    dVdt = vcat(dVdt, dVdt[end])  # add a last point same as the first
    d = tau .* dVdt .+ data
    return dVdt, -d

end

function measure_noise(data::Array{Float64,1}; threshold = 2.0, iterations = 2)
    ## Determine the base level of noise

    if iterations > 1
        med = median(data)
        stdev = std(data)
        thresh = stdev * threshold
        datam = data[(data.>med).|(data.<-med)]
        return measureNoise(datam, threshold, iterations - 1)
    else
        return std(data)
    end
end

"""
Scatter of pairs of datasets, with marginal distributions

"""

function compare_plots(x, y; labx = "", laby = "", binsx = 25, binsy = 25, ms = 2)
    println("max y: ", maximum(y))
    p2 = plot(
        x,
        y,
        seriestype = :scatter,
        markercolor = "black",
        markerstokewidth = 0.1,
        markersize = 2,
        markeralpha = 0.5,
        xlim = [0, maximum(x)],
        xlabel = labx,
        ylim = [0, maximum(y)],
        ylabel = laby,
        # subplot = 1,
        framestyle = :box,
    )
    p3 = plot(
        x,
        seriestype = :histogram,
        bins = binsx,
        orientation = :v,
        framestyle = :semi,
        xlim = [0, maximum(x)],
        ylabel = labx,
    )
    p4 = plot(
        y,
        seriestype = :histogram,
        bins = binsy,
        orientation = :h,
        framestyle = :semi,
        ylim = [0, maximum(y)],
        xlabel = laby,
    )

    return (scatter = p2, mx = p3, my = p4)
end

"""
Plot some detected measures against each other, along
with marginal histograms
events should be tuple: 
events = (dur=ev_dur*1e3, amp=ev_amp*1e12, sum=ev_sum, indx=ev_index, ipeaks=ev_peak_index)

"""
function plot_measures(x, data, crosses, events; labx = "", laby = "")
    ncrosses = crosses.ncrosses
    pcrosses = crosses.pcrosses

    default(fillcolor = :lightgrey, markercolor = :white, grid = false, legend = false)
    l = @layout [aa{0.1h}; a{0.25h}; [b{0.75w,0.15h}; c{1.1w,0.8h} d{1h,0.15w}]]
    # plot(layout = l, link = :both, size = (500, 500), margin = -10Plots.px)
    p1 = plot(x, data, linewidth = 0.5, color = "black")
    p1 = plot!(
        p1,
        x[ncrosses],
        data[ncrosses],
        seriestype = :scatter,
        markercolor = "blue",
        markerstrokewidth = 0,
        markersize = 4,
    )
    p1 = plot!(
        p1,
        x[pcrosses],
        data[ncrosses],
        seriestype = :scatter,
        markercolor = "red",
        markerstrokewidth = 0,
        markersize = 4,
    )
    p1 = plot!(
        p1,
        x[events.ipeak],
        data[events.ipeak],
        seriestype = :scatter,
        markercolor = "cyan",
        markerstrokewidth = 0,
        markersize = 4,
    )
    title = plot(
        title = @sprintf("%s Test", mode),
        grid = false,
        showaxis = false,
        yticks = false,
        xticks = false,
        bottom_margin = -50Plots.px,
    )
    hp = compare_plots(events.durs, events.amps, labx = "Dur (ms)", laby = "Amp (pA)")
    plot(title, p1, hp.mx, hp.hist, hp.my, layout = l) # , layout=grid(3, 1, heights = [0.1, 0.45, 0.35]))  # note use of show=true here - critical!
    plot!(size = (600, 600), show = true)
    show()
end

struct Event
    indices::Tuple{Int64,Int64}
    trace::Int64
    eventno::Int64
    latency::Float64
    amplitude::Float64
    onsettime::Float64
    peaktime::Float64
    risetime::Float64
    duration::Float64
    falltime::Float64
    class::String
    newclass::String
end

struct Events
    label::String
    events::Vector{Event}
end

function unpack_events(d)
    # unpack d from evts and return parallel arrays
    n = size(d.events)[1]
    trace = Vector{Int64}(undef, n)
    eventno = Vector{Int64}(undef, n)
    amp = Vector{Float64}(undef, n)
    peakt = Vector{Float64}(undef, n)
    onsettime = Vector{Float64}(undef, n)
    lat = Vector{Float64}(undef, n)  # relative to nearest previous 
    dur = Vector{Float64}(undef, n)
    rt = Vector{Float64}(undef, n)
    ft = Vector{Float64}(undef, n)
    rfratio = Vector{Float64}(undef, n)
    class = Vector{String}(undef, n)
    newclass = Vector{String}(undef, n)
    for i = 1:n
        trace[i] = d.events[i].trace
        eventno[i] = d.events[i].eventno
        amp[i] = d.events[i].amplitude
        peakt[i] = d.events[i].peaktime
        onsettime[i] = d.events[i].onsettime
        lat[i] = d.events[i].latency
        dur[i] = d.events[i].duration
        rt[i] = d.events[i].risetime
        rfratio[i] = rt[i] / d.events[i].falltime
        class[i] = d.events[i].class
        newclass[i] = d.events[i].newclass
        
    end
    return n, trace, eventno, amp, onsettime, peakt, lat, dur, rt, ft, rfratio, class, newclass
end

function events_to_dataframe(d)
    # unpack d and reformat as a DataFrame
    n, trace, eventno, amp, onset_t, peakt, lat, dur, rt, ft, rfratio, class, newclass =
        unpack_events(d)
    df = DataFrame(
        trace = trace,
        eventno = eventno,
        amp = amp,
        onsettime = onset_t,
        peaktime = peakt,
        lat = lat,
        dur = dur,
        rt = rt,
        rfratio = rfratio,
        class = class,
        newclass = newclass
    )
    return n, df
end 

function label_events(
    tdat,
    idat,
    npks,
    ev,
    pks,
    ev_end,
    template,
    thr::Float64,
    sign::Int64,
    classifier;
    data_info = nothing,
)
    if data_info == nothing
        return nothing
    end
    evts = Events("All the Events", Vector{Event}())
    dt_seconds = mean(diff(tdat[:, 1]))

    # get the light flash times from the data_info
    stim_lats = Acq4Reader.get_stim_times(data_info)
    # println("Stimlats: ", stim_lats)
    for i = 1:size(idat)[2]
        if npks[i] == 0
            continue
        else
            for j = 1:npks[i]
                index = (j, i)
                amp = idat[pks[i][j], i] .* sign .* 1e12
                rt = (pks[i][j] .- ev[i][j]) .* dt_seconds .* 1e3
                dur = (ev_end[i][j] .- ev[i][j]) .* dt_seconds .* 1e3
                ft = (ev_end[i][j] .- pks[i][j]) .* dt_seconds .* 1e3
                evlat = tdat[ev[i][j]] * 1e3  # absolute latency of this event
                peakt = tdat[pks[i][j], i] .* 1e3
                lat = -100.0  # push latency way out of our range so later ML classifier doesn't confuse with local events

                class = "spontaneous"  # assume everyting is spontaneous
                # check to see if event falls into an evoked window
                for k = 1:size(stim_lats)[1]  # check after each stimulus
                    if (evlat .> stim_lats[k]) .&
                       (evlat .<= (stim_lats[k] + classifier.maxEvokedLatency))
                        lat = evlat - stim_lats[k]
                        # now sort direct from evoked
                        if  (lat >= classifier.minDirectLatency) & # longer than shortest direct latency
                            (lat < classifier.minEvokedLatency) &  # but shorter than potential event
                            (rt >= classifier.minDirectRisetime) &  # with a slower rise time
                            (dur >= classifier.minDirectDuration) # and a longer duration
                            # println("Set direct: ", lat, " ", dur)
                            class = "direct"  # relabel
                        elseif (lat >= classifier.minEvokedLatency) & # latency and short duration mark putative evoked
                            (lat < classifier.maxEvokedLatency) &  # logically redundantwith first if statement 
                            (dur < classifier.minDirectDuration) & # not as long as a direct event
                            (dur >= classifier.minEvokedDuration) & # but longer than noise event
                            (amp >= classifier.minEvokedAmplitude)
                            # println("set evoked: ", i, " ", j, " ", lat)
                            class = "evoked"  # assign putative class
                        end
                    end
                end

                if (class == "spontaneous") &  # see if event qualifies as sponeaneous 
                    (amp < classifier.minSpontaneousAmplitude) |
                    (dur < classifier.minSpontaneousDuration) # not direct or evoked
                    # so see if is valid "event" to put in the spontaneous class
                    # println("Assigned to spont: ", i, " ", j, " ", amp, " ", dur)
                    class = "noise"
                end
                # println("lat: ", lat, " amp: ", amp, " dur: ", dur, " class: ", class)
                e = Event(index, i, j, lat, amp, evlat, peakt, rt, dur, ft, class, class)  # note "newclass" is same as class here
                push!(evts.events, e)
            end
        end
    end
    return evts
end

"""
The data structure EventTrace holds information about
an individual event - it's type, the waveform,
the trace number it cme from, and the index for the
start and end of the event (from the zero-crossing algorithm)
Typically, this would be used as an element in a vector holding
a group of events. The vector here is loaded up by extract_events,
which gets all of the events for a given class into an array

"""
struct EventTrace
    class::String
    tdat::Vector{}
    idat::Vector{}
    trace::Int64
    onset::Int64
    pkend::Int64
end

"""
    extract_events)     tdat,
    idat,
    s,
    c,
    npks,
    ev,
    pks,
    ev_end,
    template,
    thr::Float64,
    sign::Int64,
    classifier;
    data_info = nothing,)
    
    Return all extracted events meeting the classifiation, with their time bases
    
"""

function extract_events(
    tdat,
    idat,
    npks,
    df,
    sign,
    classifier,
)

    dt_seconds = mean(diff(tdat[:, 1]))
    tmax = 1.0 # seconds
    extracted = Vector{EventTrace}()
    ntraces = size(idat)[2]
    for i = 1:ntraces
        eventdf = df[in([i]).(df.trace), :]
        evk = filter(r -> any(occursin.([classifier], r.class)), eventdf)
        # evk = eventdf[in("evoked").(eventdf.class), :]
        # changed_df = filter(r -> r.class != r.newclass, eventdf)
        # println("Events: ", size(evk))
        # nchev = size(changed_df)[1]  # get number of rows with changes
        pkt = evk[(evk.peaktime.<tmax*1e3), :peaktime]
        onset = evk[(evk.peaktime.<tmax*1e3), :onsettime]
        amp = sign .* evk[(evk.peaktime.<tmax*1e3), :amp]
        dur = evk[(evk.peaktime.<tmax*1e3), :dur]
        for k = 1:size(onset)[1]
            onset_i = Int64(floor(1e-3 .* onset[k] ./ dt_seconds))
            pkend_i = Int64(floor(1e-3 .* (onset[k] .+ dur[k]) ./ dt_seconds))
            evx = EventTrace(classifier,
                    tdat[onset_i:pkend_i, i].-tdat[onset_i, i],
                    idat[onset_i:pkend_i, i],
                    i,  # save trace of origin
                    onset_i, # and onset position
                    pkend_i, # and end position
                ) 
            push!(extracted, evx)
        end
    end
    return extracted
end


"""
locate events of any shape in a signal. Works by finding regions of the signal
that deviate from noise, using the area, duration and minimum peak beneath the deviation as the detection criteria.

Makes the following assumptions about the signal:
  - noise is gaussian
  - baseline is centered at 0 (high-pass filtering may be required to achieve this).
  - no 0 crossings within an event due to noise (low-pass filtering may be required to achieve this)
  - Events last more than minDuration time
  Return an array of events where each row is (start, length, sum, peak)


"""
function zero_crossings(
    data::Array{Float64,1},
    pars::NamedTuple{
        (
            :sign,
            :dt,
            :minDuration,
            :minPeak,
            :minCharge,
            :noiseThreshold,
            :checkplot,
            :extra,
        ),
        Tuple{Int64,Float64,Float64,Float64,Float64,Float64,Bool,Bool},
    },
)

    minDuration = pars.minDuration
    minPeak = pars.minPeak
    minCharge = pars.minCharge
    noiseThreshold = pars.noiseThreshold
    data = data .* convert(Float64, pars.sign)  # flip sign if needed
    ncrosses = findall((data[1:end-1] .> 0) .& (data[2:end] .<= 0))  ## get indices of points crossing 0 in the right direction
    if sign == -1
        ncrosses = findall((data[1:end-1] .> 0) .& (data[2:end] .<= 0))  ## get indices of points crossing 0 in the right direction
        ncrosses .+= 1 # adjust array indexing
        pcrosses = findall((data[1:end-1] .<= 0) .& (data[2:end] .> 0))  ## get indices of points crossing 0 in the right direction
    else
        pcrosses = findall((data[1:end-1] .> 0) .& (data[2:end] .<= 0))  ## get indices of points crossing 0 in the right direction
        pcrosses .+= 1 # adjust array indexing
        ncrosses = findall((data[1:end-1] .<= 0) .& (data[2:end] .> 0))  ## get indices of points crossing 0 in right direction
    end
    x = collect(range(0.0, size(data)[1] - 1, step = 1.0))

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
    ev_dur = zeros(Float64, size(ncrosses)[1])
    k = 0
    l_n = size(ncrosses)[1]
    p_n = size(pcrosses)[1]
    # println(l_n, " ", p_n)
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
    if sign == -1
        n = l_n
    else
        n = p_n
    end
    for i = 1:n
        if np
            dwin = data[ncrosses[i]:pcrosses[i]]  # aligned
        else
            if i > (p_n - 1)
                # println("I exceeds p_n")
                continue  # skip as no end in sight
            end
            dwin = data[ncrosses[i]:pcrosses[i+1]]  # always go from n to p
        end
        if (minDuration > 0) & (size(dwin)[1] * pars.dt < minDuration)
            # println("rejected on length")
            continue
        end
        # println(minimum(dwin), " ", pars.minPeak)
        if (pars.minPeak > 0.0) & (maximum(dwin) < pars.minPeak)
            # println("rejected on min peak")
            continue
        end
        if (minCharge > 0.0) & (sum(dwin) * size(dwin)[1] * pars.dt < minCharge)
            # println("rejected on sum ", minCharge, " ", sum(dwin))
            continue
        end
        if (noiseThreshold > 0.0) & (maximum(dwin) < noise_value)
            continue
        end
        # should deal with noise threshold here
        k += 1
        ev_index[k] = ncrosses[i]
        ev_amp[k], pki = findmax(dwin)
        ev_peak_index[k] = convert(Int64, pki) + ncrosses[i]
        ev_sum[k] = sum(dwin)
        # ev_peak_index[k] += ev_index[k] - 1
        ev_dur[k] = size(dwin)[1] * pars.dt
    end
    ev_index = ev_index[1:k]
    ev_amp = ev_amp[1:k]
    ev_sum = ev_sum[1:k]
    ev_peak_index = ev_peak_index[1:k]
    ev_dur = ev_dur[1:k]

    crosses = (ncrosses = ncrosses, pcrosses = pcrosses)
    events = (
        durs = ev_dur .* 1e3,
        amps = ev_amp .* 1e12,
        sums = ev_sum,
        indx = ev_index,
        ipeak = ev_peak_index,
    )
    if pars.checkplot
        plot_measures(x, data, crosses, events, labx = "Dur (ms)", laby = "Amp (ms)")
    end
    crit = zeros(Float64, size(data)[1])
    crit[ev_index] .= 1.0
    if pars.extra  # return other information
        return data, crit, events
    else
        return data, crit
    end

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

function filter_data(tdat, idat; lpf = 3000, hpf = 0)
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

function make_test_waveforms(;
    N = 1,
    rate = 5.0,
    maxt = 5.0,
    lpf = 3000.0,
    taus = [1e-3, 3e-3],
    sign = -1,
)
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
        Random.seed!(42 + i)
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
        w_psc[:, i] = filter_data(t_psc, wv, lpf = lpf)
    end
    return t_psc, w_psc, dt_seconds, template

end


function identify_events(wave, crit, thr, dt, sign; thresh = 3.0, ev_win = 5e-3)
    #=
    find events crossing threshold - these are onsets
    Operates on a single trace
    =#
    # println("thr, sign: ", thr, " ", sign, " ", maximum(crit))
    ev = findall((crit[1:end-1] .> thr))
    if ev == []  # no events found!
        return [], [], thr
    end
    ev = insert!(ev, 1, 0)
    evn = findall(diff(ev) .> 1) .+ 1
    ev = ev[evn]
    swin = floor(Int, ev_win / dt)

    # find points for maximum and end of event
    pks = zeros(Int64, (size(ev)[1]))
    ev_end = zeros(Int64, (size(ev)[1]))
    maxn = size(ev)[1]
    k = 1
    for iev in ev
        if iev + swin < size(wave)[1]
            if sign == -1
                iend = findfirst((wave[iev:end-1] .<= 0) .& (wave[iev+1:end] .> 0)) # define end as first crossing
                ipk = argmin(wave[iev:(iev+iend)])
            else
                iend = findfirst((wave[iev:end-1] .>= 0) .& (wave[iev+1:end] .< 0))
                ipk = argmax(wave[iev:(iev+iend)])
            end
            pks[k] = ipk + iev - 1
            if iend != nothing
                evt = iev + ipk + iend - 2
                ev_end[k] = minimum((size(wave)[1], evt))
            end
            k += 1
        else
            maxn -= 1  # event goes beyond end, so delete it
        end
    end
    return ev[1:maxn], pks[1:maxn], ev_end[1:maxn]
end

"""
detect_events(mode, idat, dt_seconds; parallel=false,  N=1, thresh=3.0, tau1=1, tau2=3, sign=1, zcpars = ())

Detect events using different methods set by mode (CB, AJ, ZC, RS) in an array of current traces acquired
with time intervals in dt_seconds.
If parallel is true, we will run in parallel mode; otherwise each trace is analyzed serially
thresh is the threshold level for event detection measured as the sd of the criteria array
returned from the detection routines.
tau1 and tau2 are for the template (used in CB, AJ and RS)
sign is the sign of the signal to detect (+1 or -1)
zcpars is a tuple of the parameters for the zero-crossing method if it us used.
Must be of format: zcpars=(minDuration=0.0, minPeak=5e-12, minCharge=0.0, noiseThreshold=4.)
"""
function detect_events(
    mode,
    idat,
    dt_seconds;
    parallel = false,
    thresh = 3.0,
    tau1 = 1e-3,
    tau2 = 3e-3,
    sign = 1,
    zcpars = (),
)
    n_traces = size(idat)[2]
    nwave = size(idat)[1]
    taus = [tau1, tau2]
    tmax = 5 * taus[2]
    t_template, template =
        make_template(taus = taus, tmax = tmax, dt_seconds = dt_seconds, sign = sign)

    # println("make template: ", size(template))
    # println(taus)
    # println(tmax)
    # println(dt_seconds)
    ntemplate = nwave - size(template)[1]

    if mode == "CB"
        method = clements_bekkers
        pars = (t = template, dt = dt_seconds)
    elseif mode == "AJ"
        method = aj_deconvolve
        pars = (template = template, tau = [tau1, tau2], dt = dt_seconds)
    elseif mode == "RS"
        method = rs_deconvolve
        pars = (tau = tau1, dt = dt_seconds)  # tau instead of template
    elseif mode == "ZC"
        method = zero_crossings
        if zcpars == ()  # make sure that zcpars is populated
            zcpars = (
                minDuration = 0.0,
                minPeak = 5e-12,
                minCharge = 0.0,
                noiseThreshold = 4.0,
                checkplot = true,
                extra = false,
            )
        end
        pars = (
            sign = sign,
            dt = dt_seconds,
            minDuration = zcpars.minDuration,
            minPeak = zcpars.minPeak,
            minCharge = zcpars.minCharge,
            noiseThreshold = zcpars.noiseThreshold,
            checkplot = zcpars.checkplot,
            extra = zcpars.extra,
        )
        thresh = 0.5
    else
        println("Method not recognized: ", mode)
    end
    if parallel
        println(WHITE_FG, "    Event detection, parallel")
        s = SharedArray{Float64,2}((nwave, n_traces))
        c = SharedArray{Float64,2}((nwave, n_traces))
        @time @threads for i = 1:n_traces
            s[:, i], c[:, i] = method(idat[:, i], pars)
        end

    else
        println(WHITE_FG, "    Event detection, serial")
        s = Array{Float64,2}(undef, (nwave, n_traces))
        c = Array{Float64,2}(undef, (nwave, n_traces))
        @time for i = 1:n_traces
            s[:, i], c[:, i] = method(idat[:, i], pars)
        end
    end
    println(WHITE_FG, "    Detection complete") # draw the threshold line
    thr = std(c) * thresh  # all traces go into the threshold
    # println("thr: ", thr, " ", thresh)

    if parallel
        println(WHITE_FG, "    Event identification, parallel")
        evx = SharedArray{Int64,2}((nwave, n_traces))  # event onsets
        pksx = SharedArray{Int64,2}((nwave, n_traces))  # event peaks
        ev_endx = SharedArray{Int64,2}((nwave, n_traces))  # event peaks
        nppks = SharedArray{Int64,1}((n_traces))
        @time @threads for i = 1:n_traces
            eva, pksa, evenda = identify_events(idat[:, i], c[:, i], thr, dt_seconds, sign)
            nppks[i] = size(eva)[1]
            if nppks[i] > 0
                evx[1:nppks[i], i] = eva
                pksx[1:nppks[i], i] = pksa
                ev_endx[1:nppks[i], i] = evenda
            end
        end
        ev = Array{Array{Int64}}(undef, (n_traces))  # event onsets
        pks = Array{Array{Int64}}(undef, (n_traces))  # event peaks
        ev_end = Array{Array{Int64}}(undef, (n_traces))  # event peaks
        npks = Array{Int64}(undef, (n_traces))
        for i = 1:n_traces
            npks[i] = nppks[i]
            if npks[i] > 0
                ev[i] = evx[1:npks[i], i]
                pks[i] = pksx[1:npks[i], i]
                ev_end[i] = ev_endx[1:npks[i], i]
            end
        end

    else
        println(WHITE_FG, "    Event identification, serial")
        ev = Array{Array{Int64}}(undef, (n_traces))  # event onsets
        pks = Array{Array{Int64}}(undef, (n_traces))  # event peaks
        ev_end = Array{Array{Int64}}(undef, (n_traces))  # event peaks
        npks = Array{Int64,1}(undef, (n_traces))
        @time for i = 1:n_traces
            ev[i], pks[i], ev_end[i] =
                identify_events(idat[:, i], c[:, i], thr, dt_seconds, sign)
            npks[i] - size(ev[i])[1]
        end
    end
    println(WHITE_FG, "    Identification Complete")
    return s, c, npks, ev, pks, ev_end, thr

end

function plot_events(t_psc, idat, s, c, npks, ev, pks, ev_end, template, thr, sign)
    nwave = size(idat)[1]
    N_tests = size(idat)[2]

    thrline = [thr, thr]
    thrx = [0, maximum(t_psc[1:nwave, 1])]
    p1 = ""
    p2 = ""
    p3 = ""
    p4 = ""
    p5 = ""
    dt_seconds = mean(diff(t_psc[:, 1]))
    evamps = Vector{Float64}()
    evdurs = Vector{Float64}()
    evlats = Vector{Float64}()
    for i = 1:N_tests
        if npks[i] == 0
            continue
        end
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
                p1,
                t_psc[:, i],
                idat[:, i],
                #color = "black",
                palette = :Dark2_5,
                linewidth = 0.5,
                # title = "EPSC",
                # titlefontsize = 10,
                # legend = false,
            )
        end
        # plot onsets of events
        p1 = plot!(
            p1,
            t_psc[ev[i], i],
            idat[ev[i], i],
            seriestype = :scatter,
            markercolor = "yellow",
            markerstrokecolor = "black",
            markerstrokewidth = 0.1,
            markersize = 2.5,
        )
        # plot peaks of events
        p1 = plot!(
            p1,
            t_psc[pks[i], i],
            idat[pks[i], i],
            seriestype = :scatter,
            markercolor = "red",
            markerstrokecolor = "red",
            markerstrokewidth = 0.1,
            markersize = 2.5,
        )
        # plot end of events
        p1 = plot!(
            p1,
            t_psc[ev_end[i], i],
            idat[ev_end[i], i],
            seriestype = :scatter,
            markercolor = "blue",
            markerstrokecolor = "blue",
            markerstrokewidth = 0.1,
            markersize = 2.5,
        )

        # plot scale value returned
        if i == 1
            p2 = plot(
                t_psc[1:size(s)[1], i],
                s[:, i],
                linecolor = "cyan",
                linewidth = 0.3,
                title = "Scale",
                titlefontsize = 10,
                legend = false,
            )
        else
            p2 = plot!(
                p2,
                t_psc[1:size(s)[1], i],
                s[:, i],
                linecolor = "cyan",
                linewidth = 0.3,
                title = "Scale",
                titlefontsize = 10,
                legend = false,
            )
        end

        # plot criteria value (used for thresholding) returned
        if i == 1
            p3 = plot(
                t_psc[1:size(c)[1], i],
                c[:, i],
                #linecolor = "red",
                linewidth = 0.3,
                palette = :Dark2_5,
                title = "Criteria",
                titlefontsize = 10,
                legend = false,
            )
            p3 = plot!(p3, thrx, thrline, linecolor = "green", legend = false)
        else
            p3 = plot!(
                p3,
                t_psc[1:size(c)[1], i],
                c[:, i],
                #linecolor = "red",
                linewidth = 0.3,
                palette = :Dark2_5,
                title = "Criteria",
                titlefontsize = 10,
                legend = false,
            )
            p3 = plot!(p3, thrx, thrline, linecolor = "green", legend = false)
        end
        if i == 1
            if template != nothing
                p4 = plot(
                    t_psc[1:size(template)[1]],
                    template,
                    linecolor = "blue",
                    linewidth = 0.3,
                    title = "Template",
                    titlefontsize = 10,
                    legend = false,
                )
            end
        else
            if template != nothing
                p4 = plot!(
                    p4,
                    t_psc[1:size(template)[1]],
                    template,
                    linecolor = "blue",
                    linewidth = 0.3,
                    title = "Template",
                    titlefontsize = 10,
                    legend = false,
                )
            end
        end
        append!(evamps, idat[pks[i], i] .* sign .* 1e12)
        append!(evdurs, (ev_end[i] .- ev[i]) .* dt_seconds .* 1e3)
        evl = t_psc[ev[i], i]
        # append!(evlats, t_psc[ev[i], i])
        # if i == 1
        #     p5 = plot(
        #         t_psc[1:size(template)[1]],
        #         template,
        #         linecolor = "blue",
        #         linewidth = 0.3,
        #         title = "Amp vs. Width",
        #         titlefontsize = 10,
        #         legend = false,
        #     )
        # else
        #     p5 = plot!(p5,
        #         t_psc[1:size(template)[1]],
        #         template,
        #         linecolor = "blue",
        #         linewidth = 0.3,
        #         title = "Amp vs. Width",
        #         titlefontsize = 10,
        #         legend = false,
        #     )
        # end
    end

    hp = compare_plots(evdurs, evamps; labx = "Durs", laby = "Amps", binsx = 25, binsy = 25)

    title = plot(
        title = @sprintf("%s Test", mode),
        grid = false,
        showaxis = false,
        yticks = false,
        xticks = false,
        bottom_margin = -50Plots.px,
    )
    if template != nothing
        l = @layout(
            [
                a{0.1h}
                b{0.25h}
                c{0.1h}
                d{0.1h}
                e [f1{0.75w,0.15h}; f2{1.0w,0.8h} f3{1h,0.15w}]
            ]
        )
        plot(
            title,
            p1,
            p2,
            p3,
            p4,
            hp.mx,
            hp.scatter,
            hp.my,
            layout = l, # grid(5, 1, heights = [0.1, 0.25, 0.25, 0.25, 0.15, 0.15]),
        ) #, 0.30, 0.30]))

    else
        l = @layout(
            [
                a{0.1h}
                b{0.25h}
                c{0.1h}
                d{0.1h}
                [f1{0.75w,0.15h}; f2{1.0w,0.8h} f3{1h,0.15w}]
            ]
        )
        plot(
            title,
            p1,
            p2,
            p3,
            hp.mx,
            hp.scatter,
            hp.my,
            layout = l, # grid(5, 1, heights = [0.1, 0.25, 0.25, 0.25, 0.15, 0.15]),
        ) #, 0.30, 0.30]))
    end

    u = plot!(size = (600, 600))
    return u
    # Plots.savefig("CB_test.pdf")
end


function testdetector(
    mode = "CB";
    parallel = false,
    N = 1,
    thresh = 3.0,
    tau1 = 1 * 1e-3,
    tau2 = 3 * 1e-3,
    sign = 1,
    zcpars = (),
)
    if zcpars == ()
        zcpars = (
            minDuration = 0.0,
            minPeak = 5e-12,
            minCharge = 0.0,
            noiseThreshold = 4.0,
            checkplot = false,
            extra = false,
        )
    end
    N_tests = N
    t_psc, idat, dt_seconds, template =
        make_test_waveforms(N = N_tests, taus = [tau1, tau2], sign = sign)
    s, c, npks, ev, pks, ev_end, thr = detect_events(
        mode,
        idat,
        dt_seconds,
        parallel = parallel,
        thresh = thresh,
        tau1 = tau1,
        tau2 = tau2,
        sign = sign,
        zcpars = zcpars,
    )
    # u = plot_events(t_psc, idat, s, c, npks, ev, pks, ev_end, template, thr, sign)
    u = LSPSPlotting.plot_event_distributions(
        t_psc,
        idat,
        s,
        c,
        npks,
        ev,
        pks,
        ev_end,
        template,
        thr,
        sign,
    )
    # savefig("mini_events.pdf")
end

end
