module LSPSAnalysis
using LsqFit
using Statistics
using Printf

using DSP
using Base.Threads
# ENV["MPLBACKEND"] = "MacOSX"
using Plots
pyplot()


include("Acq4Reader.jl")
include("MiniAnalysis.jl")

export LSPS_read_and_plot

function ZScore(tdat, idat; baseline = [0, 0.1], score_win = [0.8, 0.85])
    ntraces = size(tdat)[2]
    bpts = findall((tdat[:, 1] .>= baseline[1]) .& (tdat[:, 1] .< baseline[2]))
    spts = findall((tdat[:, 1] .>= score_win[1]) .& (tdat[:, 1] .< score_win[2]))
    mpost = mean(idat[spts, :], dims = 1)
    mpre = mean(idat[bpts, :], dims = 1)
    zs = abs.((mpost .- mpre) ./ std(idat[bpts, :], dims = 1))
    return zs
end

function LSPSAnalyzer(tdat, vdat, idat; timewindow = [0.1, 0.6])
    x = 1

end

function Sample_and_Hold(tdat, idat; twin = [0.055, 0.056])
    #=
    Remove rapid transient artifacts from traces
    The data in twin in each trace in idat is replaced by the
    single data value in the point immediately preceding the window
        
    =#

    ipts = findall((tdat[:, 1] .>= twin[1]) .& (tdat[:, 1] .< twin[2]))
    hold_pt = ipts[1] - 1
    for i = 1:size(idat)[2]
        idat[ipts, i] .= idat[hold_pt, i]
    end
    return idat
end

# Signal processing functions for arrays (rather than SeisData or SeisChannel)
# taked from SeisNoise.jl

"""
    detrend!(X::AbstractArray{<:AbstractFloat})
Remove linear trend from array `X` using least-squares regression.
"""
function detrend!(X::AbstractArray{<:AbstractFloat})
    T = eltype(X)
    N = size(X,1)

    # create linear trend matrix
    A = similar(X,T,N,2)
    A[:,2] .= T(1)
    A[:,1] .= range(T(0),T(1),length=N)
    # create linear trend matrix
    R = transpose(A) * A

    # do the matrix inverse for 2x2 matrix
    # this is really slow on GPU
    Rinv = inv(Array(R)) |> typeof(R)
    factor = Rinv * transpose(A)

    # remove trend
    X .-= A * (factor * X)
    return nothing
end
detrend(A::AbstractArray{<:AbstractFloat}) = (U = deepcopy(A);
        detrend!(U);return U)
# detrend!(R::RawData) = detrend!(R.x)
# detrend(R::RawData) = (U = deepcopy(R); detrend!(U.x); return U)
# detrend!(C::CorrData) = detrend!(C.corr)
# detrend(C::CorrData) = (U = deepcopy(C); detrend!(U.corr); return U)


"""
    demean!(A::AbstractArray{<:AbstractFloat})
Remove mean from array `A`.
"""
function demean!(A::AbstractArray{<:AbstractFloat}; dims=1)
      μ = mean(A,dims=dims)
      A .-= μ
  return nothing
end
demean(A::AbstractArray{<:AbstractFloat}; dims=1) = (U = deepcopy(A);
       demean!(U,dims=dims);return U)
# demean!(R::RawData) = demean!(R.x)
# demean(R::RawData) = (U = deepcopy(R); demean!(U.x); return U)
# demean!(C::CorrData) = demean!(C.corr)
# demean(C::CorrData) = (U = deepcopy(C); demean!(U.corr); return U)

# ===========================


function Filter_Data(tdat, idat; lpf = 3000, hpf = 0)
    #=
    Digitally filter the data in dat with either a
    lowpass filter (lpf) or an hpf filter (hpf)
    Values are in Hz if tdat is in seconds.
    Returns filtered idat
    =#
    fs = 1 / mean(diff(tdat[:, 1]))  # get fs in Hz
    println("Fs: ", fs)
    if hpf == 0
        responsetype = Lowpass(lpf; fs = fs)
    else
        responsetype = Bandpass(hpf, lpf; fs = fs)
    end
    designmethod = Butterworth(4)
    idat = filt(digitalfilter(responsetype, designmethod), idat)
    return idat
end

function LSPS_plot_one_trace(i, p_I, tdat, idat, npks, peaks, offset, tmax, ylims, lc, lw)
    #=
    Plot one trace, including the markers on the peaks of identified events.
    
    i : trace number to plot
    p_I : the plotting object (created when i = 1, and returned)
    tdat : 1D array of time
    idat : 1D array of current
    pks : 1D array of indices to the peaks of detected events
    offset : 1D array (trials) of offsets for the traces
    tmax, top_lims : plotting parameters
    lc, lw : line color and line width for this trace (set by
    the calling routine based on some criteria, such as ZScore).
    =#
    

    

    if i == 1
        tmax = maximum(tdat)
        p_I = plot(
            tdat,
            idat .- offset[i],
            linewidth = lw,
            color = lc,
            xlim = [0., tmax],# xlabel="Dur (ms)",
            ylim = ylims, #ylabel="Amp (pA)", 
            # subplot = 1,
            framestyle = :none,
            legend=false,
        )

    else
        p_I = plot!(
            p_I,
            tdat,
            idat .- offset[i],
            linewidth = lw,
            color = lc,
            legend=false,
        )
    end
    if npks > 0
        p_I = plot!(
        p_I,
        tdat[peaks],
        idat[peaks] .- offset[i],
        seriestype = :scatter,
        markercolor = :red,
        markerstrokecolor = :red,
        markerstrokewidth = 0.5,
        markersize = 2,
        legend=false,
    )
    end
    return p_I
end

function LSPS_StackPlot(tdat, idat, data_info, npks, events, peaks, above_zthr; mode="Undef", saveflag = false)
    #=
    Make a stacked set of plots
    tdat : 2D array of (time x trial)
    idat : 2D array of (time x trial)
    data_info : the data info from the file/protocol
    top_lims : ylimits on the axes
    above_zthr : array of boolean based on AScore > some value
    saveflag : set true to write file to disk
    =#
    
    vspc = 50e-12
    twin = [0.0, 1.0]
    tmax = maximum(twin)
    ipts = findall((tdat[:, 1] .>= twin[1]) .& (tdat[:, 1] .< twin[2]))
    maxpts = maximum(ipts)

    # @time evts, pks, thr = MiniAnalysis._events(
    #     tdat[ipts, :],
    #     idat[ipts, :],
    #     taus = [4e-3, 15e-3],
    #     thresh = 3.0,
    # )
    # println("Event finding done")
    # println("Events: ", evts)
    # println("thr: ", thr)

    # i = 1
    # p_I = plot(tdat[ipts,1], idat[ipts,1] .- (1*vspc), xlims=(0, tmax), ylims=top_lims, legend=false, w=0.5, linecolor="black")
    # p_I = plot!(p_I, tdat[pks[1][:], i], idat[pks[1][:], i] .- (1*vspc),
    #     seriestype=:scatter, markercolor="red", markersize=1.5)
    ntraces = size(tdat)[2]
    p_I = 0
    offset  = Array{Float64, 1}(undef, (ntraces))
    for i = 1:ntraces
        offset[i] = i*vspc
    end
    top_lims = maximum(offset)
    # print("toplims: ", top_lims)
    bot_lims = minimum(offset)
    ylilms = [bot_lims, top_lims]
    # establish first trace

    i = 1
    lc = "black"
    lw = 0.3
    if i in above_zthr
        lc = "red"
        lw = 0.7
    end

    pks = peaks[i][findall(peaks[i] .< maxpts)]
    p_I = ""
    p_I = LSPS_plot_one_trace(i, p_I, tdat[ipts,i], idat[ipts,i], npks[i], pks, offset, tmax, ylims, lc, lw)
    
    # now do the rest of the traces
    @timed @threads for i = 2:ntraces
        lc = "black"
        lw = 0.3
        if i in above_zthr
            lc = "red"
            lw = 0.7
        end
        if npks[i] > 0
            pks = peaks[i][findall(peaks[i] .< maxpts)]
        else
            pks = []
        end
        p_I = LSPS_plot_one_trace(i, p_I, tdat[ipts,i], idat[ipts,i], npks[i], pks, offset, tmax, ylims, lc, lw)
    end
    avg = mean(idat[ipts, above_zthr], dims = 2)
    rawavg = mean(idat[ipts, :], dims = 2)
    # p_avg = plot(tdat[ipts, 1], rawavg * 1e12, w = 0.2, linecolor = "gray")
    # p_avg = plot!(p_avg, tdat[ipts, 1], avg * 1e12, w = 0.5, linecolor = "blue")
    title = plot(
        title = @sprintf("%s Test", mode),
        grid = false,
        showaxis = false,
        yticks = false,
        xticks = false,
        bottom_margin = -50Plots.px,
    )
    l = @layout([a{0.1h}; b])
    PX = plot(
        title,
        p_I,
        # layout = l,
        layout = l, # grid(5, 1, heights = [0.1, 0.25, 0.25, 0.25, 0.15, 0.15]),
    ) #, 0.30, 0.30]))
    plot!(PX, size = (600, 800))
    # if saveflag
    #     save("LSPS_test.pdf")
    # else
    #     show()
    # end
    return PX
end

function LSPS_read_and_plot(filename; fits = true, saveflag = false, mode="AJ")
    #=
    Read an HDF5 file, do a little analysis on the data
        -- ivs and fitting
    and plot the result
    =#

    tdat, idat, vdat, data_info = Acq4Reader.read_hdf5(filename)
    # top_lims, bot_lims = Acq4Reader.get_lims(data_info["clampstate"]["mode"])
    maxt = 1.0
    dt_seconds = 1.0 / data_info["DAQ.Primary"]["rate"]
    @printf(
        "Data lengths: Time=%d  I=%d  V=%d  [# traces = %d]  Duration: %8.3f sec\n",
        size(tdat)[1],
        size(idat)[1],
        size(vdat)[1],
        size(tdat)[2],
        maximum(tdat[:, 1])
    )
    ntr =  size(tdat)[2] # 20
    imax = Int64(maxt / dt_seconds)
    println("imax: ", imax)

    idat = idat[1:imax, 1:ntr]
    tdat = tdat[1:imax, 1:ntr]
    
    demean!(idat)
    detrend!(idat)
    
    idat = Sample_and_Hold(tdat, idat)
    idat = Filter_Data(tdat, idat, lpf = 2500.0)
    # calculate zscores
    sign = -1
    zs = ZScore(tdat, idat.*sign, baseline = [0, 0.1], score_win = [0.9, 0.93])
    println("Max Z Score: ", maximum(zs), " of ", size(zs))
    above_zthr = findall(zs[1, :] .> 1.96)
    println("# traces above z threshod: ", size(above_zthr), " ", above_zthr)

    zcpars = (minDuration=5e-3, minPeak=5e-12, minCharge=0e-12, noiseThreshold=4.0, checkplot=false, extra=false)
    # s, c, ev, pks, thr = MiniAnalysis.detect_events(mode, idat[:,1:ntr], dt_seconds,
    #         parallel=false, thresh=3.0, tau1=1*1e-3, tau2=3*1e-3, sign=-1,
    #         zcpars = zcpars)
    s, c, npks, ev, pks, ev_end, thr = MiniAnalysis.detect_events(mode, idat, dt_seconds,
        parallel = true, thresh=3, tau1=1e-3, tau2=3e-3, sign=sign, zcpars=zcpars)
    template=nothing
    u = MiniAnalysis.plot_event_distributions(tdat, idat, s, c, npks, ev, pks, ev_end, template, thr, sign, data_info=data_info)
    savefig("mini_event_distributions.pdf")
    # PX = LSPS_StackPlot(tdat, idat, data_info, npks, ev, pks, above_zthr, mode=mode, saveflag = saveflag)
    # savefig("test.pdf")
end

end
