module LSPSAnalysis
using LsqFit
using Statistics
using Printf

using DSP
using Base.Threads

include("Acq4Reader.jl")
include("MiniAnalysis.jl")
include("EPSC_lda.jl")
include("LSPSPlotting.jl") 

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



function LSPS_read_and_plot(filename; fits = true, saveflag = false, mode="AJ")
    #=
    Read an HDF5 file, do a little analysis on the data
        -- ivs and fitting
    and plot the result
    =#
    println("reading file")
    @time tdat, idat, vdat, data_info = Acq4Reader.read_hdf5(filename)
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
    ntr =   size(tdat)[2] # 20
    imax = Int64(maxt / dt_seconds)
    println("imax: ", imax)

    idat = idat[1:imax, 1:ntr]
    tdat = tdat[1:imax, 1:ntr]
    
    println("Preparing data")
    demean!(idat)
    detrend!(idat)
    idat = Sample_and_Hold(tdat, idat)
    idat = Filter_Data(tdat, idat, lpf = 2500.0)
    # calculate zscores
    sign = -1
    println("computing ZScores")
    @time zs = ZScore(tdat, idat.*sign, baseline = [0, 0.1], score_win = [0.9, 0.93])
    println("    Max Z Score: ", maximum(zs), " of ", size(zs))
    above_zthr = findall(zs[1, :] .> 1.96)
    println("#     traces above z threshod: ", size(above_zthr), " ", above_zthr)

    zcpars = (minDuration=5e-3, minPeak=5e-12, minCharge=0e-12, noiseThreshold=4.0, checkplot=false, extra=false)
    # s, c, ev, pks, thr = MiniAnalysis.detect_events(mode, idat[:,1:ntr], dt_seconds,
    #         parallel=false, thresh=3.0, tau1=1*1e-3, tau2=3*1e-3, sign=-1,
    #         zcpars = zcpars)
    println("detecting events")
    @time s, c, npks, ev, pks, ev_end, thr = MiniAnalysis.detect_events(mode, idat, dt_seconds,
        parallel = true, thresh=3, tau1=1e-3, tau2=3e-3, sign=sign, zcpars=zcpars)
    template=nothing
    println("labeling events")
    @time events = MiniAnalysis.label_events(tdat, idat, s, c, npks, ev, pks, ev_end, template, thr, sign, data_info=data_info)
    print("No. events: ", size(events.events)[1])
    println("Plotting distributions")
    n, df = MiniAnalysis.events_to_dataframe(events)  # we lose trace info here, but file for LDA etc.
    @time u = LSPSPlotting.plot_event_distributions(df)  # distributions don't care about trace
    # savefig(u, "mini_event_distributions.pdf")
    saveflag = false
    PX = LSPSPlotting.stack_plot(tdat, idat, data_info, sign, events, above_zthr, mode=mode)
    # savefig(PX, "test.pdf")
    return u
end

end
