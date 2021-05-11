module IVAnalysis
using LsqFit
using Statistics
using Printf

ENV["MPLBACKEND"] = "MacOSX"
using Plots
using PlotlyBase
Plots.plotly()


include("Acq4Reader.jl")
include("SpikeAnalysis.jl")

export IV_read_and_plot
#=
This module provides functions for computing current-voltage relationships
from step protocols, and for fitting selected traces to a single exponential

=#
function compute_iv(tdat, vdat, idat; ss_win = [0.5, 0.6], pk_win = [0.1, 0.2])
    #=
    Compute the current-voltage relations of a data set over
    both the steady-state window and the peak wincow
        
    tdat, vdat and idate are the time, voltage and current traces 
        (dimensions are npoints, ntraces)
    ss_win is the steady-state window
    pk_win is the window for measuring the peak voltage deflection
    Returns V, I and a minimum value
    =#
    pts = findall((tdat[:, 1] .>= ss_win[1]) .& (tdat[:, 1] .< ss_win[2]))
    pts2 = findall((tdat[:, 1] .>= pk_win[1]) .& (tdat[:, 1] .< pk_win[2]))
    vm = mean(vdat[pts, :], dims = 1)
    im = mean(idat[pts, :], dims = 1)
    imp = minimum(idat[pts2, :], dims = 1)
    return vm, im, imp
end


function fit_iv(
    tdat,
    vdat,
    idat;
    ilim = [-10e-9, -0.05e-9],
    iwin = [0.1, 0.2],
    p00 = [-0.06, -0.01, 200.0],
)
    #=
    Fit multiple traces in a current-voltage data set to measure time constants

    tdat, vdat and idate are the time, voltage and current traces 
        (dimensions are npoints, ntraces)
    ilim : the min and max currents for fitting
    p00 : Initial values for the curve fit (DC, amplitude, 1/tau)
    iwin : the window over which the data will be fit
        
    =#
    # pts =  findall((tdat[:,1] .>= window[1]) .& (tdat[:,1] .< window[2]))
    ipts = findall((tdat[:, 1] .>= iwin[1]) .& (tdat[:, 1] .< iwin[2]))
    imn = mean(idat[ipts, :], dims = 1)
    imnp = findall((imn .>= ilim[1]) .& (imn .<= ilim[2]))
    nfits = size(imnp)[1]
    npts = size(ipts)[1]
    @printf("# of fits: %d\n", nfits)
    tfit = Array{Float64,2}(undef, npts, nfits)
    vfit = Array{Float64,2}(undef, npts, nfits)
    params = Array{Any,1}(undef, nfits)
    expmodel(t, p) = p[1] .+ p[2] * exp.(-t ./ p[3])
    for i = 1:nfits
        td = tdat[ipts, i]
        p0 = p00
        vd = vdat[ipts, i]
        fit = curve_fit(expmodel, td .- iwin[1], vd, p0)
        params[i] = fit.param
        @printf(
            "Params: DC= %8.2f mV A = %8.2f mV  Tau = %8.3f ms\n",
            params[i][1] * 1e3,
            params[i][2] * 1e3,
            params[i][3] * 1e3
        )
        tfit[:, i] = td
        vfit[:, i] = expmodel(td .- iwin[1], params[i])
    end
    return tfit, vfit, params
end

function IV_read_and_plot(filename, fits = true, ivs = true, analyze_spikes = true)
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
    if ivs
        vm, im, imp = IVAnalysis.compute_iv(tdat, vdat, idat)
    end
    if fits
        tfit, vfit, params =
            IVAnalysis.fit_iv(tdat, vdat, idat, iwin = [0.1, 0.15], ilim = [-1e-9, 0])
        tfit2, vfit2, params2 =
            IVAnalysis.fit_iv(tdat, vdat, idat, iwin = [0.2, 0.6], ilim = [-1e-9, 0])
        # println("params: ", params)
    end
    if analyze_spikes
        SpikeAnalysis.AnalyzeSpikes(tdat, vdat, idat)
    end
    tmax = maximum(tdat[:, 1])
    if ivs
        p0 = plot(
            hcat(vm[1, :], vm[1, :]),
            hcat(im[1, :], imp[1, :]),
            line = true,
            w = 3.0,
            m = 5,
            marker = "circle",
        )
    end
    p2 = plot(
        tdat,
        idat,
        xlims = (0, tmax),
        ylims = top_lims,
        legend = false,
        w = 0.5,
        linecolor = "black",
    )
    p1 = plot(
        tdat,
        vdat,
        xlims = (0, tmax),
        ylims = bot_lims,
        legend = false,
        w = 0.3,
        linecolor = "black",
    )
    if fits
        plot!(p1, tfit, vfit, linecolor = "red")
        plot!(p1, tfit2, vfit2, linecolor = "blue")
    end# p1d = plot(td, vd, xlims=(0, tmax), ylims=top_lims, legend=false)
    # p1f = plot(tfit, vfit, )
    l = @layout [a{0.6w} b]
    leftplots = plot(p1, p2, layout = grid(2, 1, heights = [0.75, 0.25]))
    if ivs
        plot(leftplots, p0, layout = l, legend = false)
    end
    plot!(size = (800, 500))
    # plot(p1, p2, p0, layout = grid(2, 2, heights=[0.75, 0.25]), title=("V", "I"))
end

end
