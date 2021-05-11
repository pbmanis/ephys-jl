module LSPSAnalysis
using LsqFit
using Statistics
using Printf

using DSP
using Base.Threads

# ENV["MPLBACKEND"] = "MacOSX"
using GLMakie

# using PyPlot
# using PlotlyBase
# Plots.plotly()


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

function LSPS_plot_one_trace(ax, i, p_I, tdat, idat, pks, offset, tmax, top_lims, lc, lw)
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
    
    ### Original using Plots/plotly
    # if i == 1
    #     p_I = plot(
    #         tdat,
    #         idat .- offset[i],
    #         xlims = (0, tmax),
    #         ylims = top_lims,
    #         legend = false,
    #         w = lw,
    #         linecolor = lc,
    #     )
    #
    # else
    #     p_I = plot!(
    #         tdat,
    #         idat .- offset[i],
    #         xlims = (0, tmax),
    #         ylims = top_lims,
    #         legend = false,
    #         w = lw,
    #         linecolor = lc,
    #     )
    # end
    # p_I = plot!(
    #     p_I,
    #     tdat[pks],
    #     idat[pks] .- offset[i],
    #     seriestype = :scatter,
    #     markercolor = "red",
    #     markeralpah = 0.8,
    #     markerstrokecolor = "white",
    #     markerstrokewidth = 0.5,
    #     markerstrokealpha = 0.5,
    #     markersize = 1.5,
    #     legend=false,
    # )
    # return p_I
    
    #Using Makie
    

    if i == 1
        p_I = lines!(
            ax, tdat,
            idat .- offset[i],
            linewidth = lw,
            color = lc,
        )

    else
        p_I = lines!(
            ax, tdat,
            idat .- offset[i],
            linewidth = lw,
            color = lc,
        )
    end
    p_I = scatter!(
        ax, tdat[pks],
        idat[pks] .- offset[i],
        color = :red,
        strokecolor = :white,
        strokewidth = 0.5,
        markersize = 1.5,
    )
    return p_I
end

function LSPS_StackPlot(tdat, idat, data_info, top_lims, above_zthr; saveflag = false)
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


    @time evts, pks, thr = MiniAnalysis.CB_find_events(
        tdat[ipts, :],
        idat[ipts, :],
        taus = [4e-3, 15e-3],
        thresh = 3.0,
    )
    println("Event finding done")
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


    fig = Figure(backgroundcolor = RGBf0(0.98, 0.98, 0.98),
        resolution = (800, 1000))
    ax = Axis(fig[2,1])
    println("ax: ", ax)
    xlims!(ax, [0, tmax])
    println("toplims: ", top_lims)
    tl = top_lims * 1e-12
    ylims!(ax, (0.0, tl))
    
    # establish first trace

    i = 1
    lc = "black"
    lw = 0.3
    if i in above_zthr
        lc = "red"
        lw = 0.7
    end
    p_I = ""
    p_I = LSPS_plot_one_trace(ax, i, p_I, tdat[ipts,i], idat[ipts,i],  pks[i], offset, tmax, top_lims, lc, lw)
    
    # now do the rest of the traces
    @timed @threads for i = 2:ntraces
        lc = "black"
        lw = 0.3
        if i in above_zthr
            lc = "red"
            lw = 0.7
        end
        p_I = LSPS_plot_one_trace(ax, i, p_I, tdat[ipts,i], idat[ipts,i], pks[i], offset, tmax, top_lims, lc, lw)
    end
    avg = mean(idat[ipts, above_zthr], dims = 2)
    rawavg = mean(idat[ipts, :], dims = 2)
    # p_avg = plot(tdat[ipts, 1], rawavg * 1e12, w = 0.2, linecolor = "gray")
    # p_avg = plot!(p_avg, tdat[ipts, 1], avg * 1e12, w = 0.5, linecolor = "blue")
    # title = plot(
    #     title = "LSPS",
    #     grid = false,
    #     showaxis = false,
    #     yticks = false,
    #     xticks = false,
    #     bottom_margin = -50Plots.px,
    # )
    # plot(
    #     title,
    #     p_I,
    #     p_avg,
    #     layout = grid(3, 1, heights = [0.05, 0.85, 0.1]),
    #     legend = false,
    # )
    # plot!(size = (600, 800))
    # plot(p1, p2, p0, layout = grid(2, 2, heights=[0.75, 0.25]), title=("V", "I"))
    current_figure()
    # if saveflag
    #     save("LSPS_test.pdf")
    # else
    #     show()
    # end
    
end

function LSPS_read_and_plot(filename; fits = true, saveflag = false)
    #=
    Read an HDF5 file, do a little analysis on the data
        -- ivs and fitting
    and plot the result
    =#

    tdat, idat, vdat, data_info = Acq4Reader.read_hdf5(filename)
    top_lims, bot_lims = Acq4Reader.get_lims(data_info["clampstate"]["mode"])
    @printf(
        "Data lengths: Time=%d  I=%d  V=%d  [# traces = %d]  Duration: %8.3f sec\n",
        size(tdat)[1],
        size(idat)[1],
        size(vdat)[1],
        size(tdat)[2],
        maximum(tdat[:, 1])
    )

    idat = Sample_and_Hold(tdat, idat)
    idat = Filter_Data(tdat, idat, lpf = 2500.0)
    # calculate zscores
    zs = ZScore(tdat, idat, baseline = [0, 0.1], score_win = [0.9, 0.93])

    println("Max Z Score: ", maximum(zs), " of ", size(zs))
    above_zthr = findall(zs[1, :] .> 1.96)
    println("# traces above z threshod: ", size(above_zthr), " ", above_zthr)

    @time LSPS_StackPlot(tdat, idat, data_info, top_lims, above_zthr, saveflag = saveflag)

end

end
