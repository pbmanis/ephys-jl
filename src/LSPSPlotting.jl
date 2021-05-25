module LSPSPlotting
using Base
import DataFrames: DataFrame, describe, select, Not
using Statistics
using Printf
# ENV["MPLBACKEND"] = "MacOSX"
using Plots
pyplot()
include("MiniAnalysis.jl")

export plot_one_trace, stack_plot
export plot_event_distribution

function plot_event_distributions(df; response_window=25.0)
    # d is the dataframe (from miniAnalysis events_to_dataframe)
    binsx = 50
    p_amp = plot(df.amp, seriestype = :histogram, bins=binsx,
        orientation = :v, 
        framestyle = :semi,
        xlim = [0, maximum(df.amp)], 
        xlabel="Amplitude (pA)",
        ylabel = "# obs",
        )
    p_dur = plot(df.dur, seriestype = :histogram, bins=binsx,
        orientation = :v, 
        framestyle = :semi,
        xlim = [0, maximum(df.dur)], xlabel="Durations (ms)",
        ylabel = "# obs",
        )
    p_rt = plot(df.rt, seriestype = :histogram, bins=binsx,
        orientation = :v, 
        framestyle = :semi,
        xlim = [0, maximum(df.rt)], xlabel="Rise Times (ms)",
        ylabel = "# obs",
        )
    p_lat = plot(df.lat[df.lat .>= 0.], seriestype = :histogram, bins=binsx,
        orientation = :v, 
        framestyle = :semi,
        xlim = [0, response_window], xlabel="Latencies (ms)",
        ylabel = "# obs",
        )
    l = @layout[a b; c d] # ; e f]
    u = plot(p_amp, p_dur, p_rt, p_lat, layout=l, show=true)
    return u
end


function plot_one_trace(i, p_I, tdat, idat, eventdf, offset, sign, ylims, lc, lw)
    #=
    Plot one trace, including the markers on the peaks of identified events.
    Used by stack_plot
    i : trace number to plot
    p_I : the plotting object (created when i = 1, and returned)
    tdat : 1D array of time
    idat : 1D array of current
    events: The event structure, but with only those events in this trace selected
    offset : 1D array (trials) of offsets for the traces
    tmax, top_lims : plotting parameters
    lc, lw : line color and line width for this trace (set by
    the calling routine based on some criteria, such as ZScore).
    =#

    tmax = maximum(tdat)
    if i == 1
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
    evk = filter(r -> any(occursin.(["evoked"], r.class)), eventdf) # eventdf[in("evoked").(eventdf.class), :]
    pkt = evk[(evk.peaktime .< tmax*1e3), :peaktime]
    amp = sign .* evk[(evk.peaktime .< tmax*1e3), :amp]
    p_I = plot!(
        p_I,
        pkt .* 1e-3,
        (amp .* 1e-12) .- offset[i],
        seriestype = :scatter,
        markercolor = :green,
        markerstrokecolor = :green,
        markerstrokewidth = 0.5,
        markersize = 3,
        legend=false,
    )
    evk = filter(r -> any(occursin.(["direct"], r.class)), eventdf) # eventdf[in("evoked").(eventdf.class), :]
    pkt = evk[(evk.peaktime .< tmax*1e3), :peaktime]
    amp = sign .* evk[(evk.peaktime .< tmax*1e3), :amp]
    p_I = plot!(
        p_I,
        pkt .* 1e-3,
        (amp .* 1e-12) .- offset[i],
        seriestype = :scatter,
        markercolor = :orange,
        markerstrokecolor = :orange,
        markerstrokewidth = 0.5,
        markersize = 4,
        legend=false,
    )
    evk = filter(r -> any(occursin.(["spontaneous"], r.class)), eventdf) # eventdf[in("evoked").(eventdf.class), :]
    pkt = evk[(evk.peaktime .< tmax*1e3), :peaktime]
    amp = sign .* evk[(evk.peaktime .< tmax*1e3), :amp]
    p_I = plot!(
        p_I,
        pkt .* 1e-3,
        (amp .* 1e-12) .- offset[i],
        seriestype = :scatter,
        markercolor = :cyan,
        markerstrokecolor = :cyan,
        markerstrokewidth = 0.5,
        markersize = 3,
        legend=false,
    )
    return p_I
end

function stack_plot(tdat, idat, data_info, sign, events, above_zthr; mode="Undef", saveflag = false)
    #=
    Make a stacked set of plots
    tdat : 2D array of (time x trial)
    idat : 2D array of (time x trial)
    data_info : the data info from the file/protocol
    top_lims : ylimits on the axes
    above_zthr : array of boolean based on AScore > some value
    saveflag : set true to write file to disk
    =#
    
    println("Plotting stacs of traces")
    vspc = 50*1e-12
    twin = [0.0, 1.0] 
    tmax = maximum(twin)*1e3 # express in msec
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
    n, df = MiniAnalysis.events_to_dataframe(events)
    i = 1
    lc = "black"
    lw = 0.3
    if i in above_zthr
        lc = "red"
        lw = 0.7
    end

    p_I = ""
    p_I = plot_one_trace(i, p_I, tdat[ipts,i], idat[ipts,i], df[in([1]).(df.trace), :], offset,  sign, ylims, lc, lw)
    # now do the rest of the traces
    @timed  for i = 2:ntraces
        lc = "black"
        lw = 0.3
        if i in above_zthr
            lc = "red"
            lw = 0.7
        end
        p_I = plot_one_trace(i, p_I, tdat[ipts,i], idat[ipts,i], df[in([i]).(df.trace), :], offset, sign, ylims, lc, lw)
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
    plot!(PX, size = (600, 800), show=true)
    # if saveflag
        savefig("LSPS_test.pdf")
    # else
    show()
    # end
    return PX
end

end
