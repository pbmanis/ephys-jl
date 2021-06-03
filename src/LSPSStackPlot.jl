module LSPSStackPlot
using Base
import DataFrames: DataFrame, describe, select, Not
using Statistics
using Printf
# ENV["MPLBACKEND"] = "MacOSX"
#using GLMakie
using CairoMakie
CairoMakie.activate!(type = "pdf")
# using AbstractPlotting
# AbstractPlotting.__init__()
include("MiniAnalysis.jl")

export stack_plot


function stack_plot(
    df,
    tdat,
    idat,
    data_info,
    sign,
    events,
    above_zthr;
    mode = "Undef",
    figurename = Union{str, nothing} = nothing,
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

    # n, df = MiniAnalysis.events_to_dataframe(events)
    figure = Figure(resolution = (600, 800))
    ax1 = figure[1, 1] = Axis(figure, title = "Pre Treatment")
    idx = [vcat(tdat[ipts, i], NaN) for i in 1:ntraces]
    idy = [vcat(idat[ipts, i] .+ vertical_offset[i], NaN) for i in 1:ntraces]
    idx = reduce(vcat, idx)
    idy = reduce(vcat, idy)
    @time lines!(ax1, idx, idy, linewidth=0.35, overdraw=true, color=:black)

    # avg = mean(idat[ipts, above_zthr], dims = 2)
    # rawavg = mean(idat[ipts, :], dims = 2)
    # p_avg = plot(tdat[ipts, 1], rawavg * 1e12, w = 0.2, linecolor = "gray")
    # p_avg = plot!(p_avg, tdat[ipts, 1], avg * 1e12, w = 0.5, linecolor = "blue")

    stim_lats = 1e-3 .* MiniAnalysis.get_stim_times(data_info)
    stx = [vcat([stim_lats[i], stim_lats[i]], NaN) for i = 1:size(stim_lats)[1]]
    sty = [vcat([maximum(vertical_offset), 0.0], NaN)  for i = 1:size(stim_lats)[1]]
    stx = reduce(vcat, stx)
    sty = reduce(vcat, sty)
    lines!(ax1, stx, sty, linewidth=0.5, color=:lightblue, overdraw=true)
    # title = plot(
    #     title = @sprintf("%s Test", mode),
    #     grid = false,
    #     showaxis = false,
    #     yticks = false,
    #     xticks = false,
    #     bottom_margin = -50Plots.px,
    # )
    # l = @layout([a{0.1h}; b])
    # PX = plot(
    #     title,
    #     p_I,
    #     # layout = l,
    #     layout = l, # grid(5, 1, heights = [0.1, 0.25, 0.25, 0.25, 0.15, 0.15]),
    # ) #, 0.30, 0.30]))
    #
    # PX =  plot!(p_I, size = (600, 800), show = true)
    # if figurename != nothing
    #     savefig(scene, figurename)
    # end
    # display()
    save("stack_plot.pdf", figure)
    # GLMakie.display(figure)
    return figure
end



end
