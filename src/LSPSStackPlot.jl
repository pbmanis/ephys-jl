module LSPSStackPlot
using Base
import DataFrames: DataFrame, describe, select, Not
using Statistics
using Printf
# ENV["MPLBACKEND"] = "MacOSX"
#using GLMakie
using GLMakie
using CairoMakie

# whichMakie.activate!()  # interactive first
using AbstractPlotting
AbstractPlotting.__init__()
AbstractPlotting.inline!(false)
include("MiniAnalysis.jl")
include("Acq4Reader.jl")

export stack_plot2



function stack_plot2(
    df,
    tdat,
    idat,
    data_info,
    sign,
    events,
    above_zthr;
    mode = "Undef",
    figurename = "",
    maxtraces = 0,
    makie = "i", # for interactive; other is "p" for pdf output
)
    #=
    Make a stacked set of plots, using Makie
    Trick: the arrays are all concatenated so there is only one plot object for the big data
    tdat : 2D array of (time x trial)
    idat : 2D array of (time x trial)
    data_info : the data info from the file/protocol
    top_lims : ylimits on the axes
    above_zthr : array of boolean based on AScore > some value
    saveflag : set true to write file to disk
    =#
    println("   Doing stacked trace plot")
    if makie == "i"
        GLMakie.activate!()
        whichMakie = GLMakie
        println("Activated GLMakie")
    elseif makie == "p"
        resolution = (6000 ,6000)
        CairoMakie.activate!()   
        whichMakie = CairoMakie
        println(("activated CairoMakie"))
    else
        println("makie must be either i (interactive) or p (pdf)")
        return nothing
    end
    vspc = 50 * 1e-12
    twin = [0.0, 1.0]
    tmax = maximum(twin) * 1e3 # express in msec
    ipts = findall((tdat[:, 1] .>= twin[1]) .& (tdat[:, 1] .< twin[2]))
    maxpts = maximum(ipts)

    ntraces = size(tdat)[2]
    if (maxtraces > 0) &( maxtraces < ntraces)
        ntraces = maxtraces
    end
    p_I = 0
    vertical_offset = Array{Float64,1}(undef, (ntraces))
    for i = 1:ntraces
        vertical_offset[i] = i * vspc
    end
    top_lims = maximum(vertical_offset)
    # print("toplims: ", top_lims)
    bot_lims = minimum(vertical_offset)
    ylims = [bot_lims, top_lims]

    # n, df = MiniAnalysis.events_to_dataframe(events)
    figure = whichMakie.Figure(resolution = resolution)
    ax1 = figure[1, 1] = whichMakie.Axis(figure)
    figure[0, :] = whichMakie.Label(figure, figurename, textsize=10)
    figure[3,1] = buttongrid = GridLayout(tellwidth = false)

    # scene = Scene(figure)
    # cam = whichMakie.cam2D!(figure) # you can also get the cam from an existing scene
    # cam[:panbutton] = Mouse.left
    
    whichMakie.hidespines!(ax1, :t, :r, :l)
    whichMakie.hideydecorations!(ax1)
    whichMakie.hidexdecorations!(ax1, ticks=false, ticklabels=false)
    ax1.xticks= 0.0:0.2:1.0
    # println(minimum(vertical_offset), maximum(vertical_offset))
    # ax1.yticks=collect(range(0, maximum(vertical_offset); step=500e-12))
    # ax1.yminorticks = 5
    # ax1.yminorticksvisible = true
    whichMakie.lines!(ax1, [-0.06, -0.06, -0.01,], [50e-12, 0, 0,], linewidth=1.25, color=:black)
    whichMakie.annotations!("50 pA", position=Point2f0(-0.07, 25e-12), textsize=8, align=(:left, :center))
    whichMakie.annotations!("50 ms", position=Point2f0(-0.035, -20e-12), textsize=8, align=(:center, :top))
    
    
    idx = [vcat(tdat[ipts, i], NaN) for i = 1:ntraces]
    idy = [vcat(idat[ipts, i] .+ vertical_offset[i], NaN) for i = 1:ntraces]
    cmap = [vcat(zeros(size(idat[ipts, i])[1]), 0) for i = 1:ntraces]  # color map per point
    dt_seconds = mean(diff(tdat[:,1]))
    # class = "direct"
    clmap = Dict("direct" => 4, "evoked" => 6, "spontaneous" =>2, "noise" => 1, "artifact" => 9)  # for use with paired_10
    for (class, clevel) in clmap
        for i = 1:ntraces
            eventdf = df[in([i]).(df.trace), :]
            evk = filter(r -> any(occursin.([class], r.class)), eventdf) # eventdf[in("evoked").(eventdf.class), :]
            # changed_df = filter(r -> r.class != r.newclass, eventdf)
            # println("Changed events: ", size(changed_df))
            # nchev = size(changed_df)[1]  # get number of rows with changes
            pkt = evk[(evk.peaktime.<tmax*1e3), :peaktime]
            onset = evk[(evk.peaktime.<tmax*1e3), :onsettime]
            amp = sign .* evk[(evk.peaktime.<tmax*1e3), :amp]
            dur = evk[(evk.peaktime.<tmax*1e3), :dur]
            for k = 1:size(onset)[1]
                onset_i = Int64(floor(1e-3 .* onset[k] ./ dt_seconds))
                pkend_i = Int64(floor(1e-3 .* (onset[k] .+ dur[k]) ./ dt_seconds))
                cmap[i][onset_i:pkend_i] .= clevel
            end
        end
    end
    
    idx = reduce(vcat, idx)
    idy = reduce(vcat, idy)
    cmap = reduce(vcat, cmap)
    println("Preparing to plot")
    @time whichMakie.lines!(ax1, idx, idy, linewidth=0.25, color=cmap, colormap=:Paired_10)

    # avg = mean(idat[ipts, above_zthr], dims = 2)
    # rawavg = mean(idat[ipts, :], dims = 2)
    # p_avg = plot(tdat[ipts, 1], rawavg * 1e12, w = 0.2, linecolor = "gray")
    # p_avg = plot!(p_avg, tdat[ipts, 1], avg * 1e12, w = 0.5, linecolor = "blue")

    stim_lats = 1e-3 .* Acq4Reader.get_stim_times(data_info, device="Laser")
    stx = [vcat([stim_lats[i], stim_lats[i]], NaN) for i = 1:size(stim_lats)[1]]
    sty = [vcat([maximum(vertical_offset), 0.0], NaN)  for i = 1:size(stim_lats)[1]]
    stx = reduce(vcat, stx)
    sty = reduce(vcat, sty)
    whichMakie.lines!(ax1, stx, sty, linewidth=0.5, color=:lightblue, overdraw=true)
    sfx = (0.0, 1.0)
    sfy = (minimum(vertical_offset)-vspc, maximum(vertical_offset)+vspc)# title = plot(
    # buttons = buttongrid = whichMakie.Button("Reset")
    
    if whichMakie == "GLMakie"
        println("Building buttons for GLMakie")
        button_labels = ["Reset Axes", "Save  png", "save  pdf"]
        buttons = buttongrid[1, 1:3] = [whichMakie.Button(figure, label = button_labels[i]) for i in 1:3]

        on(buttons[1].clicks) do n
            whichMakie.xlims!(ax1, sfx)
            whichMakie.ylims!(ax1,  sfy)
        end
        on(buttons[2].clicks) do n
            GLMakie.activate!()
            save("stack_plot2.png", figure)
            println("PNG saved")
        end
        on(buttons[3].clicks) do n
            whichMakie.activate!()
            save("stack_plot2.pdf", figure)
            println("PDF saved")
        end
    end
    println("Which Makie: ", whichMakie)
    display(figure)
    
    if whichMakie == CairoMakie
        save("stack_plot2_cairo.png", figure)
    end
    # save("stack_plot.pdf", figure)  # when done
    return figure
end



end
