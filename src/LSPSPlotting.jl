module LSPSPlotting
using Base
using Formatting
using StatsBase
import DataFrames: DataFrame, describe, select, Not
using Statistics
using Printf
# ENV["MPLBACKEND"] = "MacOSX"
using Plots
# pyplot()

include("LSPSFitting.jl")
include("Acq4Reader.jl")

export plot_one_trace, stack_plot, plot_trace, finalize_plot
export plot_event_distribution

function plot_event_distributions(df; response_window = 25.0, figurename=Union{str, nothing} = nothing)
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
    if figurename != nothing
        savefig(u, figurename)
    end
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
    changed_df = filter(r -> r.class != r.newclass, eventdf)
    # println("Changed events: ", size(changed_df))
    nchev = size(changed_df)[1]  # get number of rows with changes
    pkt = evk[(evk.peaktime.<tmax*1e3), :peaktime]
    onset = evk[(evk.peaktime.<tmax*1e3), :onsettime]
    amp = sign .* evk[(evk.peaktime.<tmax*1e3), :amp]
    dur = evk[(evk.peaktime.<tmax*1e3), :dur]
    # classify changed events so we can mark them
    #=
    Choose from [:none, :auto, :circle, :rect, :star5, :diamond,
     :hexagon, :cross, :xcross, :utriangle, :dtriangle, 
    :rtriangle, :ltriangle, :pentagon, :heptagon, :octagon, :star4, :star6, :star7, :star8, :vline, :hline, :+, :x].
    =#
    for i = 1:size(changed_df)[1]
        markersize = 3
        markershape = :circle
        dfr = changed_df[i, :]
        # println("changed: ", i, " ", dfr.class)
        if (dfr.class == "evoked") & (dfr.newclass == "direct")
            markersize=8
            markershape = :dtriangle  # demoted
        elseif (dfr.class == "direct") & (dfr.newclass == "evoked")
            markersize=8
            markershape = :utriangle # promoted
        elseif (dfr.class == "direct") & (dfr.newclass == "spontaneous")
            markersize=8
            markershape = :xcross
        elseif (dfr.class == "evoked") & (dfr.newclass == "spontaneous")
                markersize=8
                markershape = :plus  # demoted
        else
            markersize=4
            markershape = :star5  # mark all others
        end
        p_I = plot!(
            p_I,
            [dfr.peaktime*1e-3],
            [(sign*dfr.amp*1e-12) .+ vertical_offset],
            seriestype = :scatter,
            markercolor = :red,
            markersize = markersize,
            markershape = markershape,
            legend=false,
            )
    end
    
    
    # decorate with a dot at the peak using the passed color for the event type
    p_I = plot!(
        p_I,
        pkt .* 1e-3,
        (amp .* 1e-12) .+ vertical_offset,
        seriestype = :scatter,
        markercolor = markercolor,
        markerstrokecolor = markercolor,
        markeralpha = 0.35,
        markerstrokewidth = 0.5,
        markersize = 2.5,
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
            idat[onset_i:pkend_i] .+ vertical_offset,
            color = linecolor,
            linewidth = 1.0 * linewidth,
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
            idat .+ vertical_offset,
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
            idat .+ vertical_offset,
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
        linecolor = :red,
        linewidth = linewidth*1.5,
        class = "evoked",
        markercolor = :red,
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
        linecolor = :cyan,
        linewidth = linewidth*0.75,
        class = "spontaneous",
        markercolor = :cyan,
    )
    # p_I = decorate_and_paint!(
    #     p_I,
    #     tdat,
    #     idat,
    #     vertical_offset,
    #     eventdf,
    #     sign,
    #     linecolor = :magenta,
    #     linewidth = linewidth,
    #     class = "noise",
    #     markercolor = :magenta,
    # )

    return p_I
end


function fit_and_plot_events(p_raw, p_sub, x, y, color; plotting::Bool=true, mindy::Float64=100.0)

    # fits = Vector{Float64}(undef, 4)
    # for i = 1:4
    #     minflag, y0, efit = LSPSFittingfit_direct(x, y, n=i, mindy=mindy)
    #     rms = rmsd(y, y0)
    #     fits[i] = rms*1e12
    #     # println("i: ", i, "  rms: ", rms)
    # end
    # bestn = argmin(fits)
    bestn = 3
    minflag, y0, efit = LSPSFitting.fit_direct(x, y, n=bestn, mindy=mindy)
    
    if !plotting
        return p_raw, p_sub, y0
    end

    if isnothing(p_raw)
        p_raw = plot(x, y, linecolor=color, legend=false)
        p_raw = plot!(p_raw, x, y0, linecolor=color, legend=false, linestyle=:dash)
    else
        p_raw = plot!(p_raw, x, y, linecolor=color, legend=false)
        p_raw = plot!(p_raw, x, y0, linecolor=color, legend=false, linestyle=:dash)
        # println("was not nothing")
    end
    if isnothing(p_sub)
        p_sub = plot(x, y .- y0, linecolor=color, legend=false)
    else
        p_sub = plot!(p_sub, x, y .- y0, linecolor=color, legend=false)
    end
    return p_raw, p_sub, y0
end


function finalize_fitted_plot(p1, p2)
    if (p1 == nothing) | (p2 == nothing)
        return nothing
    end
    l = @layout [a; b] # ; e f]
    xlabel!(p1, "T (sec)")
    ylabel!(p1, "I (A)")
    xlabel!(p2, "T (sec)")
    ylabel!(p2, "I (A)")
    PX =  plot(p1, p2, layout=l, size = (600, 600), show = true)
    return PX
end


function stack_plot(
    df,
    tdat,
    idat,
    data_info,
    sign,
    events,
    above_zthr;
    mode = "Undef",
    figtitle = "",
    figurename = Union{str, nothing} = nothing,
    maxtraces = 0,
    makie = "" # ignored - just for compatability with LSPSStackPlot.stack_plot2
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
    println("   Doing stacked trace plot")
    vspc = 50.0 * 1e-12
    twin = [0.0, 1.0]
    tmax = maximum(twin) * 1e3 # express in msec
    ipts = findall((tdat[:, 1] .>= twin[1]) .& (tdat[:, 1] .< twin[2]))
    maxpts = maximum(ipts)
    ntraces = size(tdat)[2]
    if (maxtraces > 0) & (maxtraces < ntraces)
        ntraces = maxtraces
    end
    vertical_offset = Array{Float64,1}(undef, (ntraces))
    for i = 1:ntraces
        vertical_offset[i] = (i-1) * vspc
    end
    top_lims = vspc * (ntraces+1)
    bot_lims = - vspc
    ylims = [bot_lims, top_lims]
    println("ylims: ", ylims)
    
    p_I = nothing
    @timed for i = 1:ntraces
        # println("Plotting trace: ", i)
        lc = :black
        lw = 0.25
        if i in above_zthr
            lc = :red
            lw = 0.5
        end
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

    stim_lats = 1e-3 .* Acq4Reader.get_stim_times(data_info, device="Laser")
    for i = 1:size(stim_lats)[1]
        p_I = plot!(
            p_I,
            [(stim_lats[i], stim_lats[i])],
            [(ylims[1], ylims[2])],
            linewidth = 0.5,
            linecolor = :blue,
        )
    end
    p_I = plot!(p_I, ylims = ylims)
    labels = Dict("Evoked" => :red, "Direct" => :orange, "Spontaneous" => :gray, "Noise" => "magenta")
    i = 0
    for (class, color) in labels
        p_I = plot!(p_I, [0.05+0.1*(i-1)], [ylims[1]-1.1*vspc], linewidth=2, label=class, color=color)
        i += 1
    end
    title = plot(
        title = @sprintf("%s :: %s", mode, figurename),
        grid = false,
        showaxis = false,
        yticks = false,
        xticks = false,
        bottom_margin = -50Plots.px,
    )
    l = @layout([a{0.1h}; b])
    p_I = plot(
        title, 
        p_I,
        # layout = l,
        layout = l, # grid(5, 1, heights = [0.1, 0.25, 0.25, 0.25, 0.15, 0.15]),
    ) #, 0.30, 0.30]))

    PX =  plot!(p_I, size = (600, 800), show = true)
    if figurename != nothing
        savefig(PX, figurename)
    end
    return PX
end

end
