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

function plot_event_distributions(df; response_window = 25.0)
    # d is the dataframe (from miniAnalysis events_to_dataframe)
    binsx = 50
    p_amp = plot(
        df.amp,
        seriestype = :histogram,
        bins = binsx,
        orientation = :v,
        framestyle = :semi,
        xlim = [0, maximum(df.amp)],
        xlabel = "Amplitude (pA)",
        ylabel = "# obs",
    )
    p_dur = plot(
        df.dur,
        seriestype = :histogram,
        bins = binsx,
        orientation = :v,
        framestyle = :semi,
        xlim = [0, maximum(df.dur)],
        xlabel = "Durations (ms)",
        ylabel = "# obs",
    )
    p_rt = plot(
        df.rt,
        seriestype = :histogram,
        bins = binsx,
        orientation = :v,
        framestyle = :semi,
        xlim = [0, maximum(df.rt)],
        xlabel = "Rise Times (ms)",
        ylabel = "# obs",
    )
    p_lat = plot(
        df.lat[df.lat.>=0.0],
        seriestype = :histogram,
        bins = binsx,
        orientation = :v,
        framestyle = :semi,
        xlim = [0, response_window],
        xlabel = "Latencies (ms)",
        ylabel = "# obs",
    )
    l = @layout [a b; c d] # ; e f]
    u = plot(p_amp, p_dur, p_rt, p_lat, layout = l, show = true)
    return u
end

function decorate_and_paint!(
    p_I,
    tdat::Vector,
    idat::Vector,
    vertical_offset::Float64,
    eventdf::DataFrame,
    sign::Int64;
    linecolor::Union{Symbol,String} = :green,
    markercolor::Union{Symbol,String} = :green,
    linewidth::Float64 = 1.0,
    class::String = "evoked",
)
    #=
    Decorate traces with a dot to indicate the detected peak, 
    then paint the selected part of the trace with a color to indicate
    the event onset/duration and classification
        
    =#
    tmax = maximum(tdat)
    dt_seconds = mean(diff(tdat))
    evk = filter(r -> any(occursin.([class], r.class)), eventdf) # eventdf[in("evoked").(eventdf.class), :]
    pkt = evk[(evk.peaktime.<tmax*1e3), :peaktime]
    onset = evk[(evk.peaktime.<tmax*1e3), :onsettime]
    amp = sign .* evk[(evk.peaktime.<tmax*1e3), :amp]
    dur = evk[(evk.peaktime.<tmax*1e3), :dur]
    # decorate with a dot
    p_I = plot!(
        p_I,
        pkt .* 1e-3,
        (amp .* 1e-12) .- vertical_offset,
        seriestype = :scatter,
        markercolor = markercolor,
        markerstrokecolor = markercolor,
        markerstrokewidth = 0.5,
        markersize = 3,
        legend = false,
    )

    # paint events according to their classification
    for k = 1:size(onset)[1]
        onset_i = Int64(floor(1e-3 .* onset[k] ./ dt_seconds))
        pkend_i = Int64(floor(1e-3 .* (onset[k] .+ dur[k]) ./ dt_seconds))
        # println("onset: ", onset_i, " pkend: ", pkend_i)
        p_I = plot!(
            p_I,
            tdat[onset_i:pkend_i],
            idat[onset_i:pkend_i] .- vertical_offset,
            color = linecolor,
            linewidth = 2.5 * linewidth,
            legend = false,
        )
    end
    return p_I
end


function plot_one_trace(
    p_I,
    tdat,
    idat;
    vertical_offset::Float64=0.0,
    eventdf::DataFrame,
    sign::Int64=1,
    linecolor::Union{Symbol,String} = :black,
    linewidth::Float64 = 0.5,
)
    #=
    Plot one trace, including the markers on the peaks of identified events.
    Used by stack_plot
    p_I : the plotting object (created when i = 1, and returned)
    tdat : 1D array of time
    idat : 1D array of current
    offset : 1D array (trials) of offsets for the traces
    events: The event structure, but with only those events in this trace selected
    tmax, top_lims : plotting parameters
    lc, lw : line color and line width for this trace (set by
    the calling routine based on some criteria, such as ZScore).
    =#

    tmax = maximum(tdat)
    dt_seconds = mean(diff(tdat))
    # plot traces simply
    if isnothing(p_I)
        p_I = plot(
            tdat,
            idat .- vertical_offset,
            linewidth = linewidth,
            color = linecolor,
            # xlim = [0.0, tmax],# xlabel="Dur (ms)",
            # ylim = ylims, #ylabel="Amp (pA)",
            # subplot = 1,
            framestyle = :none,
            legend = false,
        )
    else
        p_I = plot!(
            p_I,
            tdat,
            idat .- vertical_offset,
            linewidth = linewidth,
            color = linecolor,
            legend = false,
        )


    end

    p_I = decorate_and_paint!(
        p_I,
        tdat,
        idat,
        vertical_offset,
        eventdf,
        sign,
        linecolor = :green,
        linewidth = linewidth,
        class = "evoked",
        markercolor = :green,
    )
    p_I = decorate_and_paint!(
        p_I,
        tdat,
        idat,
        vertical_offset,
        eventdf,
        sign,
        linecolor = :orange,
        linewidth = linewidth,
        class = "direct",
        markercolor = :orange,
    )
    p_I = decorate_and_paint!(
        p_I,
        tdat,
        idat,
        vertical_offset,
        eventdf,
        sign,
        linecolor = :yellow,
        linewidth = linewidth,
        class = "spontaneous",
        markercolor = :yellow,
    )

    return p_I
end

function plot_trace(p_I, x, y)
    if isnothing(p_I)
        p_I = plot(x, y, legend=false)
        println("isnothing")
    else
        p_I = plot!(p_I, x, y, legend=false)
        println("was not nothing")
    end
    return p_I
end

function stack_plot(
    tdat,
    idat,
    data_info,
    sign,
    events,
    above_zthr;
    mode = "Undef",
    saveflag = false,
)
    #=
    Make a stacked set of plots
    tdat : 2D array of (time x trial)
    idat : 2D array of (time x trial)
    data_info : the data info from the file/protocol
    top_lims : ylimits on the axes
    above_zthr : array of boolean based on AScore > some value
    saveflag : set true to write file to disk
    =#

    println("Plotting stacks of traces")
    vspc = 50 * 1e-12
    twin = [0.0, 1.0]
    tmax = maximum(twin) * 1e3 # express in msec
    ipts = findall((tdat[:, 1] .>= twin[1]) .& (tdat[:, 1] .< twin[2]))
    maxpts = maximum(ipts)

    ntraces = size(tdat)[2]
    p_I = 0
    vertical_offset = Array{Float64,1}(undef, (ntraces))
    for i = 1:ntraces
        vertical_offset[i] = i * vspc
    end
    top_lims = maximum(vertical_offset)
    # print("toplims: ", top_lims)
    bot_lims = minimum(vertical_offset)
    ylims = [bot_lims, top_lims]

    n, df = MiniAnalysis.events_to_dataframe(events)
    println("ntraces: ", ntraces)
    p_I = nothing
    @timed for i = 1:ntraces
        lc = :black
        lw = 0.3
        if i in above_zthr
            lc = :red
            lw = 0.7
        end
        # p_I = plot_trace(p_I, tdat[ipts, i], idat[ipts, i] .+ vertical_offset[i])
        p_I = plot_one_trace(
            p_I,
            tdat[ipts, i],
            idat[ipts, i],
            vertical_offset = vertical_offset[i],
            eventdf = df[in([i]).(df.trace), :],
            sign=sign,
            linecolor = lc,
            linewidth = lw,
        )

    end

    # avg = mean(idat[ipts, above_zthr], dims = 2)
    # rawavg = mean(idat[ipts, :], dims = 2)
    # p_avg = plot(tdat[ipts, 1], rawavg * 1e12, w = 0.2, linecolor = "gray")
    # p_avg = plot!(p_avg, tdat[ipts, 1], avg * 1e12, w = 0.5, linecolor = "blue")

    stim_lats = 1e-3 .* MiniAnalysis.get_stim_times(data_info)
    for i = 1:size(stim_lats)[1]
        p_I = plot!(
            p_I,
            [(stim_lats[i], stim_lats[i])],
            [(0.0, -maximum(vertical_offset))],
            linewidth = 0.5,
            linecolor = :skyblue,
        )
    end

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

    PX =  plot!(p_I, size = (600, 800), show = true)
    return PX
end

end
