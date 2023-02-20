__precompile__()
module IVAnalysis
using LsqFit
using Statistics
using Printf
using ArraysOfArrays
using ForwardDiff
using ArgParse


ENV["MPLBACKEND"] = "MacOSX"
using Plots

using Plots.Measures
# inspectdr()
# using PlotlyBase
# Plots.plotly()
gr()
using FilePathsBase

include("Acq4Reader.jl")
include("SpikeAnalysis.jl")

export IV_read_and_plot
#=
This module provides functions for computing current-voltage relationships
from step protocols, and for fitting selected traces to a single exponential

=#
# function parse_commandline()
# 	p = ArgParseSettings()
# 	@add_arg_table p begin
# 		"filename"
# 			help = "Positional: filename"
# 			required = true
# 		"--nospikes", "--n"
# 			help = "skip spike analysis"
# 			action = :store_true
# 	end
# 	return parse_args(p)
# end
#
# function main()
# 	parsed_args = parse_commandline()
# 	println("Parsed Args")
# 	for (arg.val) in parsed_args
# 		println("  $arg -> $val")
# 	end
# end

function compute_iv(tdat, vdat, idat; ss_win=[0.5, 0.6], pk_win=[0.1, 0.2])
    #=
    Compute the current-voltage relations of a data set over
    both the steady-state window and the peak wincow

    tdat, vdat and idate are the time, voltage and current traces 
        (dimensions are npoints, ntraces)
    ss_win is the steady-state window (times in seconds)
    pk_win is the window for measuring the peak voltage deflection (times in seconds)
    Returns V, I and a minimum value
    =#
    pts = findall((tdat[:, 1] .>= ss_win[1]) .& (tdat[:, 1] .< ss_win[2]))
    pts2 = findall((tdat[:, 1] .>= pk_win[1]) .& (tdat[:, 1] .< pk_win[2]))
    vm = mean(vdat[pts, :], dims=1)
    im = mean(idat[pts, :], dims=1)
    imp = minimum(idat[pts2, :], dims=1)
    return vm, im, imp
end


function fit_iv(
    tdat,
    vdat,
    idat;
    ilim=[-10e-9, -0.05e-9],
    iwin=[0.1, 0.2],
    p00=[-0.06, -0.01, 0.005],
    lb=[-0.200, -0.05, 0.0005],
    ub=[0.01, 0.100, 0.5],
)
    #=
    Fit multiple traces in a current-voltage data set to measure time constants

    tdat, vdat and idate are the time, voltage and current traces 
        (dimensions are npoints, ntraces, times are in seconds)
    ilim : the min and max currents for fitting (in Amperes)
    iwin : the window over which the data will be fit (times in seconds)
    p00 : Initial values for the curve fit (DC, amplitude, 1/tau)

    =#

    # pts =  findall((tdat[:,1] .>= window[1]) .& (tdat[:,1] .< window[2]))
    ipts = findall((tdat[:, 1] .>= iwin[1]) .& (tdat[:, 1] .< iwin[2]))
    imn = mean(idat[ipts, :], dims=1)
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
        fit = curve_fit(expmodel, td .- iwin[1], vd, p0, lower=lb, upper=ub)
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
    return tfit, vfit, imn, params
end

function get_spikes(tdat, vdat, idat, sppars; iwin=[0.1, 1.0], ilim=[0.0, 10e-9])
    #=
    get times and counts

    tdat, vdat and idate are the time, voltage and current traces 
        (dimensions are npoints, ntraces, times in seconds, currents in Amperes,
        voltages in Volts)
    sppars :  spike parameters, passed to SpikeAnalysis.AnalyzeSpikes:
        Use SpikeDetectParameters function to fill this structure.
    iwin : the window over which spikes will be accepted (times in seconds)
    ilim : the min and max currents for counting spikes (current in Amperes)
    =#

    # pts =  findall((tdat[:,1] .>= window[1]) .& (tdat[:,1] .< window[2]))
    ntraces = size(tdat)[2]
    ipts = findall((tdat[:, 1] .>= iwin[1]) .& (tdat[:, 1] .< iwin[2]))
    ih = mean(idat[1:10, :], dims=1)
    imn = mean(idat[ipts, :], dims=1) .- ih  # incorporate holding current...  
    imnp = findall((imn .>= ilim[1]) .& (imn .<= ilim[2]))
    imnpa = [a[2] for a in imnp]
    nspktraces = size(imnp)[1]
    npts = size(ipts)[1]
    @printf("# of eligible traces: %d\n", nspktraces)
    # Array{Array{Int64}}(undef, (n_traces))
    spk_times = Array{Array{Float64}}(undef, (ntraces))
    spk_vpeak = Array{Array{Float64}}(undef, (ntraces))
    spk_indices = Array{Array{Int64}}(undef, (ntraces))
    spk_inj = Array{Float64,1}(undef, (ntraces))
    for i = 1:ntraces
        if i in imnpa
            itrace = i  # get from cartesian indexint
            # println("Itrace: ", itrace)
            # println(size(tdat), size(vdat), size(idat))
            spk_times[i], spk_vpeak[i], spk_indices[i] = SpikeAnalysis.AnalyzeSpikes(
                tdat[:, itrace],
                vdat[:, itrace],
                idat[:, itrace],
                pars=sppars,
            )
            spk_inj[i] = imn[i]
        end

    end
    # for i = 1:ntraces
    # 	if isassigned(spk_times, i)
    # 		println("i = ", i, "  nspk: ", size(spk_times[i]))
    # 	end
    # end
    return spk_times, spk_vpeak, spk_indices, spk_inj
end

function IV_read_and_plot(filename, fits=true, ivs=true, analyze_spikes=true)
    #=
    Read an HDF5 file (acq4 file), do a little analysis on the data
        -- ivs and fitting
    and plot the result
    =#

    spike_vstep = 0.085 # V steps between traces for I > 0

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
    # println(data_info["DAQ.Command"])
    # println(data_info)
    pulse_times = data_info["MultiClamp1.pulse_pars"]
    pstart = pulse_times[1]
    pdur = pulse_times[1] + pulse_times[2]
    ss_win = [(pulse_times[2] - 0.1 * pulse_times[2]) + pstart, pdur]
    pk_win = [pstart, pstart + 0.1 * pulse_times[2]]
    println("ss_win: ", ss_win, " ", "pkwin: ", pk_win)
    #println(pulse_times)
    if ivs
        vm, im, imp = IVAnalysis.compute_iv(tdat, vdat, idat, ss_win=ss_win, pk_win=pk_win)
    end
    if fits
        tfit, vfit, imean, params =
            IVAnalysis.fit_iv(tdat, vdat, idat, iwin=[pulse_times[1], pulse_times[1] + 0.1], ilim=[-1e-9, 0])
        tfit2, vfit2, imean, params2 =
            IVAnalysis.fit_iv(tdat, vdat, idat, iwin=[pulse_times[1] + 0.05, pulse_times[1]+0.5], ilim=[-1e-9, -5e-11])
        # println("params: ", params)
    end
    if analyze_spikes
        sppars = SpikeAnalysis.SpikeDetectParameters()
        spk_times, spk_vpeak, spk_indices, spk_inj = IVAnalysis.get_spikes(
            tdat,
            vdat,
            idat,
            sppars,
            iwin=[pulse_times[1], pulse_times[1] + pulse_times[2]],
            ilim=[0.0, 10e-9],
        )
        spk_count = Array{Int64}(undef, (size(spk_inj)))
        for i = 1:size(spk_inj)[1]
            if isassigned(spk_indices, i)
                spk_count[i] = size(spk_indices[i])[1]
            else
                spk_count[i] = 0
            end
        end
    end
    voffset = Array{Float64}(undef, (size(vdat)[2]))
    lasti = 1
    tspikes = Vector{Float64}()
    vspikes = Vector{Float64}()

    first_spiketrace = false
    for i = 1:size(vdat)[2]
        idat0 = mean(idat[1:10, i])
        if (imean[i] <= idat0) | (spk_count[i] == 0)
            voffset[i] = 0.0
            lasti = i+1
        else
            voffset[i] = (i - lasti) * spike_vstep
            vdat[:, i] = vdat[:, i] .+ voffset[i]

            if isassigned(spk_vpeak, i) # add markers for spikes
                append!(vspikes, spk_vpeak[i] .+ voffset[i])
                append!(tspikes, spk_times[i])
            end
        end
    end


    #
    # Data is ready to plot... 

    tmax = maximum(tdat[:, 1])
    # plot the IV relationship
    if ivs
        p_ssIV = plot(
            hcat(im[1, :], imp[1, :]) * 1e9,
            hcat(vm[1, :], vm[1, :]) * 1e3,
            line=true,
            w=3.0,
            m=3.5,
            marker="circle",
        )
    end
    # Plot the FI relationship
    if analyze_spikes
        p_FI = plot(
            spk_inj * 1e9,
            spk_count,
            line=true,
            w=3.0,
            m=3.5,
            marker="circle"
        )
    end

    # Plot the current command traces

    p_Itraces = plot(
        tdat * 1e3,
        idat * 1e9,
        xlims=(0, tmax * 1e3),
        ylims=(top_lims * 1e9, -top_lims * 1e9),
        legend=false,
        w=0.5,
        linecolor="black",
    )

    fnsplit = split(filename, "/")
    fn = "File: " * fnsplit[end-3] * "/" * fnsplit[end-2] * "/" * fnsplit[end-1] * "/" * fnsplit[end]

    # plot the voltage traces
    ylim_min = minimum(vdat)
    ylim_max = maximum(vdat)
    p_Vtraces = plot(
        tdat .* 1e3,
        vdat .* 1e3,
        xlims=(0, tmax * 1e3),
        # ylims = (ylim_min, ylim_max),
        legend=nothing,
        w=0.2,
        linecolor="darkgrey",
        title=fn,
        titlefontsize=9,
    )
    # if fits
    #     p_Fits_taum = plot(
    #         tfit .* 1e3,
    #         vfit .* 1e3,
    #         linecolor="red",
    #         w=0.5,
    #     )
    # end
    p3plot = false
    p_ISI = 0
    if analyze_spikes
        plot!(
            p_Vtraces,  # put on voltage traces
            tspikes .* 1e3,
            vspikes .* 1e3,
            seriestype=:scatter,
            w=0,
            m=2.0,
            marker="circle",
            markercolor="red",
            markerstrokewidth=0,
        )
        # plot isi vs latency
        first_plot = true
        for i = 1:size(spk_times)[1]
            if !isassigned(spk_times, i)
                continue
            end
            if length(spk_times[i]) < 2
                continue
            end
            u = length(spk_times[i])
            if first_plot == true
                p_ISI = plot(
                    spk_times[i][1:u-1] * 1e3,
                    diff(spk_times[i]) * 1e3,
                    marker="circle",
                    markerstrokewidth=0,
                    m=2.5,
                    line=true,
                )
                first_plot = false
                plot!(xlabel="Latency (ms)", ylabel="ISI (ms)", p_ISI)
            else
                plot!(
                    p_ISI,
                    spk_times[i][1:u-1] * 1e3,
                    diff(spk_times[i]) * 1e3,
                    marker="circle",
                    markerstrokewidth=0,
                    m=2.5,
                    line=true,
                )
            end
        end
        p3plot = true
    end
    plot!(xlabel="T (ms)", ylabel="V (mV)", p_Vtraces)
    plot!(xlabel="T (ms)", ylabel="I (nA)", p_Itraces)
    plot!(xlabel="I (nA)", ylabel="V (mV)", p_ssIV)
    plot!(xlabel="I (nA)", ylabel="Spikes", p_FI)
    if fits
        plot!(p_Vtraces, tfit.*1e3, vfit.*1e3, w=0.35, color="red", line=:solid)
        plot!(p_Vtraces, tfit2.*1e3, vfit2.*1e3, w=0.05, color="darkred")
    end# p1d = plot(td, vd, xlims=(0, tmax), ylims=top_lims, legend=false)
    # p1f = plot(tfit, vfit, )
    l = @layout [a{0.6w} b]
    leftplots = plot(p_Vtraces, p_Itraces, layout=grid(2, 1, heights=[0.9, 0.1]))
    if p3plot
        rightplots = plot(p_ssIV, p_FI, p_ISI, layout=grid(3, 1, heights=[0.33, 0.33, 0.33]))
    else
        rightplots = plot(p_ssIV, p_FI, layout=grid(2, 1, heights=[0.33, 0.33]))
    end
    plot(leftplots, rightplots, layout=l, legend=false, margin=3mm, rightmargin=8mm)
    plot!(size=(800, 800))
    # plot(p1, p2, p0, layout = grid(2, 2, heights=[0.75, 0.25]), title=("V", "I"))
end

# main()
end
