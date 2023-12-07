module LSPSAnalysis
using LsqFit
using Statistics
using Printf
using Crayons
using Crayons.Box
using DSP
using Base.Threads
using Revise
using Plots
using Gtk

ENV["PYTHON"] = "/Users/pbmanis/Desktop/Python/ephys/ephys_venv/bin/python"
include("Acq4Reader.jl")
include("MiniAnalysis.jl")
include("EPSC_lda.jl")
include("LSPSPlotting.jl")
include("LSPSStackPlot.jl")
include("LSPSFitting.jl")

export LSPS_read_and_plot, aone

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

function sample_and_hold(tdat, idat; twin = [0.055, 0.056])
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
    N = size(X, 1)

    # create linear trend matrix
    A = similar(X, T, N, 2)
    A[:, 2] .= T(1)
    A[:, 1] .= range(T(0), T(1), length = N)
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
    detrend!(U);
    return U)

"""
    demean!(A::AbstractArray{<:AbstractFloat})
Remove mean from array `A`.
"""
function demean!(A::AbstractArray{<:AbstractFloat}; dims = 1)
    μ = mean(A, dims = dims)
    A .-= μ
    return nothing
end
demean(A::AbstractArray{<:AbstractFloat}; dims = 1) = (U = deepcopy(A);
demean!(U, dims = dims);
return U)
# demean!(R::RawData) = demean!(R.x)
# demean(R::RawData) = (U = deepcopy(R); demean!(U.x); return U)
# demean!(C::CorrData) = demean!(C.corr)
# demean(C::CorrData) = (U = deepcopy(C); demean!(U.corr); return U)

# ===========================

function filter_data_notches(
    tdat,
    idat;
    notchfreqs = [60.0],
    Q = 20.0,
    bw = 5.0,
    QScale = true,
)
    #=
    Digitally filter the data in idat with a sequence of notch filters
    of specified q and sample frequency
    Values are in Hz if tdat is in seconds.
    Returns filtered idat
    =#
    fs = 1 / mean(diff(tdat[:, 1]))  # get fs in Hz
    w0 = notchfreqs ./ (fs / 2.0)
    # if QScale
    #     bw = w0[0]/Q
    #     Qf = (w0 ./ bw).^sqrt(2.0)
    # else
    #     Qf = Q .* ones(size(notchfreqs)[1])
    # end
    # println("Qf: ", Qf)
    # println("w0: ", w0)
    # println("bw: ", bw)
    # bw = ones(size(notchfreqs)[1]) .* range(5, 20, length=length(notchfreqs)[1])
    bw =
        ones(size(notchfreqs)[1]) .*
        exp10.(range(0.301, 1.301, length = length(notchfreqs)[1]))

    for i = 1:size(notchfreqs)[1]
        responsetype = Bandstop(notchfreqs[i] - bw[i], notchfreqs[i] + bw[i], fs = fs)
        designmethod = Butterworth(16)
        idat = filt(digitalfilter(responsetype, designmethod), idat)
    end
    return idat
end



function filter_data_bandpass(tdat, idat; lpf = 3000, hpf = 0)
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
    designmethod = Butterworth(4)
    idat = filt(digitalfilter(responsetype, designmethod), idat)
    return idat
end



function LSPS_analyze(;fits = true, 
    saveflag = false, mode = "AJ", maxtr=0, makie="i", 
    mindy=10.0, adjust=false)
    #=
    present user with a file selector interface; exit if canceled.
    otherwise try to do analysis on the directory
    =#
    filename = ""
    while filename == ""
        filename = open_dialog("Select Dataset Folder", action=GtkFileChooserAction.SELECT_FOLDER)
        println("File: ", filename)
        if (filename != "") & isdir(filename)
            LSPS_read_and_plot(filename=filename, fits=fits, saveflag=saveflag, mode=mode,
            maxtr = maxtr, makie=makie, mindy=mindy,
            adjust=adjust)
            filename = ""
        elseif filename == ""
            return
        end
    end
end
    
function LSPS_read_and_plot(;filename::String="", fits = true, 
    saveflag = false, mode = "AJ", maxtr=0, makie="i", 
    mindy=10.0, adjust=false)
    println("Reading Data File, V0.5")
    println("Filename: ", filename)
    
    #=
    Read an HDF5 file, do a little analysis on the data
        -- ivs and fitting
    and plot the result
    =#
    
    @time tdat, idat, vdat, data_info = Acq4Reader.read_hdf5(filename)
    # top_lims, bot_lims = Acq4Reader.get_lims(data_info["clampstate"]["mode"])
    maxt = 1.0
    dt_seconds = 1.0 / data_info["DAQ.Primary"]["rate"]
    @printf(
        "    Data lengths: Time=%d  I=%d  V=%d  [# traces = %d]  Duration: %8.3f sec\n",
        size(tdat)[1],
        size(idat)[1],
        size(vdat)[1],
        size(tdat)[2],
        maximum(tdat[:, 1])
    )

    ntr =  size(tdat)[2] # 20
    
    if (maxtr > 0) & (maxtr <= ntr)
        ntr = maxtr
    end
    imax = Int64(maxt / dt_seconds)
    if imax > size(tdat)[1]
        imax = size(tdat)[1]
        maxt = imax * dt_seconds 
    end
    println("ntr: ", ntr, " maxtr: ", maxtr, " dt_sec: ", dt_seconds, " imax: ", imax, " maxt: ", maxt)
    idat = idat[1:imax, 1:ntr]
    tdat = tdat[1:imax, 1:ntr]

    println("Preparing data (detrending, artifact suppression, filtering)")
    demean!(idat)
    detrend!(idat)
    notchfreqs = [
        60.0, #120.0, 180.0,
        #240.0,
        4000.0,
    ]
    idat = sample_and_hold(tdat, idat)
    # idat = filter_data_notches(tdat, idat, notchfreqs=notchfreqs)
    idat = filter_data_bandpass(tdat, idat, lpf = 2500.0)
    # calculate zscores
    sign = -1
    println("Computing ZScores")
    @time zs = ZScore(tdat, idat .* sign, baseline = [0, 0.1], score_win = [0.9, 0.93])
    println("    Max Z Score: ", maximum(zs), " of ", size(zs))
    zthr = 1.96
    above_zthr = findall(zs[1, :] .> zthr)
    println("    # traces above ZScore threshold of ", zthr, " : ", size(above_zthr))

    zcpars = (
        minDuration = 1e-3,
        minPeak = 5e-12,
        minCharge = 0e-12,
        noiseThreshold = 4.0,
        checkplot = false,
        extra = false,
    )
    
    classifier = (  # parameters for classifying events
        minEvokedAmplitude = 3.0,
        minEvokedDuration = 0.3,
        minEvokedCharge = 0.,  # not implemented
        minEvokedLatency = 3.0,
        maxEvokedLatency = 15.0,
        minDirectLatency = 0.0,
        minDirectDuration = 10.0,
        minDirectRisetime = 2.0,
        minSpontaneousAmplitude = 3.0,
        minSpontaneousDuration = 0.3
    )
    # s, c, ev, pks, thr = MiniAnalysis.detect_events(mode, idat[:,1:ntr], dt_seconds,
    #         parallel=false, thresh=3.0, tau1=1*1e-3, tau2=3*1e-3, sign=-1,
    #         zcpars = zcpars)
    println("Detecting events")
    @time s, c, npks, ev, pks, ev_end, thr = MiniAnalysis.detect_events(
        mode,
        idat,
        dt_seconds,
        parallel = true,
        thresh = 2,
        tau1 = 3e-3,
        tau2 = 10e-3,
        sign = sign,
        zcpars = zcpars,
    )
    template = nothing
    println("Labeling events (classifcation)")
    @time events = MiniAnalysis.label_events(
        tdat,
        idat,
        npks,
        ev,
        pks,
        ev_end,
        template,
        thr,
        sign,
        classifier,
        data_info = data_info,
    )
    println("    Number of events: ", size(events.events)[1])

    n, df = MiniAnalysis.events_to_dataframe(events)  # we lose trace info here, but file for LDA etc.
    # newdf = EPSC_LDA.epsc_lda(df)  # do classification and re-estimation

    # println("Plotting distributions")
    # @time u = LSPSPlotting.plot_event_distributions(df, figurename="mini_event_distributions.pdf")  #distributionsdon't care about trace

    println("Extracting direct events", )
    extracted = MiniAnalysis.extract_events(tdat, idat, npks, df, sign, "direct")
    println("Number of extracted events of class direct: ", size(extracted))
    nextracted = size(extracted)[1]
    p_ex = nothing
    p_diff = nothing
    println("N extracted: ", nextracted)
    if nextracted > 0
        CList = cgrad(:Paired_10, nextracted, categorical=true)
    else
        CList = []
    end
    println("mindy: ", mindy)
    for i in 1:nextracted
        ex = extracted[i]
        p_ex, p_diff, yfit = LSPSPlotting.fit_and_plot_events(p_ex, p_diff, ex.tdat, ex.idat, 
                CList[2], plotting=adjust, mindy=mindy)
        ysub = ex.idat .- yfit
        tr = ex.trace
        idat[ex.onset:ex.pkend, ex.trace] .= ysub
    end
    if adjust
        px = LSPSPlotting.finalize_fitted_plot(p_ex, p_diff)
        return
    end
    #
    # println("Plotting traces and annotations", )
    # @time PX = LSPSPlotting.stack_plot(
    #     df,
    #     tdat,
    #     idat,
    #     data_info,
    #     sign,
    #     events,
    #     above_zthr,
    #     mode = mode,
    #     figurename = "stack_plot.pdf"
    # )
    # return PX, u

    # repeat analysis, this time with AJ method now that directs are removed
    println("Detecting events, with AJ")
    mode = "AJ"
    @time s, c, npks2, ev2, pks2, ev_end2, thr2 = MiniAnalysis.detect_events(
        mode,
        idat,
        dt_seconds,
        parallel = true,
        thresh = 2.5,
        tau1 = 3e-3,
        tau2 = 10e-3,
        sign = sign,
        zcpars = zcpars,
    )
    template = nothing
    println("Labeling events (classifcation)")
    @time events2 = MiniAnalysis.label_events(
        tdat,
        idat,
        npks2,
        ev2,
        pks2,
        ev_end2,
        template,
        thr2,
        sign,
        classifier,
        data_info = data_info,
    )
    println("    Number of events: ", size(events2.events)[1])

    n, df = MiniAnalysis.events_to_dataframe(events2)  # we lose trace info here, but file for LDA etc.
    # newdf = EPSC_LDA.epsc_lda(df)  # do classification and re-estimation
    


    splitname = splitpath(filename)
    figtitle = joinpath(splitname[end-3:end]...)
    # @time PX = LSPSPlotting.stack_plot(
    @time PX = LSPSPlotting.stack_plot(
        df,
        tdat,
        idat,
        data_info,
        sign,
        events2,
        above_zthr,
        mode = mode,
        figurename = "stack_plot1.pdf",
        figtitle = figtitle,
        maxtraces = ntr,
        makie = makie,
    )
    return PX
end

end