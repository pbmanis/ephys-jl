__precompile__()
module Acq4Reader
using Statistics
using HDF5
using Plots
using Printf
using Base.Threads
using ArraysOfArrays
using ElasticArrays
using Crayons.Box
using Distributed
using SharedArrays  # important for parallel/looping
# using ProgressMeter
using InteractiveUtils

using PyCall

# pgc = pyimport("pyqtgraph.configfile")  # use python configuration reader...

# Theoretically, the following should work, but we encounter an error:
# "[NameError: name 'array' is not defined]")
# on the scanner target array structure with lots of "array([ 0.00157848,  0.003536  ])"
scriptdir = @__DIR__
pushfirst!(PyVector(pyimport("sys")."path"), scriptdir)
pgc = pyimport("python.configfile")  # use a local python routine to read the config files

export read_hdf5, get_lims, get_stim_times

pyplot()



function get_subdirs(base_dir)
    # read the subdirectories of the current base_dir
    # only return those that are protocol sweeps
    subdirs = readdir(base_dir)
    deleteat!(subdirs, subdirs .== ".index")
    deleteat!(subdirs, subdirs .== ".DS_Store")
    return subdirs
end

function read_hdf5(filename)
    #=
    Read the protocol data from the specified metaarray file,
    treating it as an HDF5 file (which it is)

    Return the data, and some "standard" y limits.
    New parallel version, 5/24/2021 pbm
    =#
    device = "MultiClamp1.ma"
    laser_device = "Laser-Blue-raw.ma"
    photodiode_device = "Photodiode.ma"

    sweeps = get_subdirs(filename)
    idat = Array{Float64}(undef)
    vdat = Array{Float64}(undef)
    tdat = Array{Float64}(undef)
    nsweeps = size(sweeps)[1]
    println(WHITE_FG, "    Nsweeps: ", nsweeps)
    # temporary allocation so we don't lose scope on arrays
    nwave = 1
    s_idat = SharedArray{Float64,2}((nwave, nsweeps))
    s_vdat = SharedArray{Float64,2}((nwave, nsweeps))
    s_tdat = SharedArray{Float64,2}((nwave, nsweeps))

    first = true
    mode = ""
    clampstate = ""
    data_info = ""
    println("    Number of threads: ", Threads.nthreads())
    # note we wet this up for threading, but that causes memory errors...
    # kept same structure here though.
    @time for s = 1:nsweeps
        sweep = sweeps[s]
        time, data, data_info = read_one_sweep(filename, sweep, device)
        if time == false
            continue
        end
        # will be VC, IC or IC=0 (may depend on age of acquisition code)
        sweep_mode = String(data_info["clampstate"]["mode"])
        if first
            mode = sweep_mode
        else
            if mode != sweep_mode
                throw(
                    ErrorException(
                        "Mode changed from ",
                        _mode,
                        " to ",
                        sweep_mode,
                        "inside protocol",
                    ),
                )
            end
        end

        indx1, indx2 = get_indices(data_info)

        if first  # we don't know allocation size until reading the first run
            nwave = size(data)[1]
            s_idat = SharedArray{Float64,2}((nwave, nsweeps))
            s_vdat = SharedArray{Float64,2}((nwave, nsweeps))
            s_tdat = SharedArray{Float64,2}((nwave, nsweeps))
            first = false
            # s_tdat[:,s] = time
            # s_idat[:,s] = data[:, indx1]
            # s_vdat[:,s] = data[:, indx2]
        end
        s_tdat[:, s] = time
        s_idat[:, s] = data[:, indx1]
        s_vdat[:, s] = data[:, indx2]
        # println("s: ", s)
    end
    if mode == ""
        finalize(s_tdat)
        finalize(s_idat)
        finalize(s_vdat)
        @everywhere GC.gc()
        println("mode: ", mode)
        throw(ErrorException("No acquisition mode found"))
    end
    idat = deepcopy(s_idat)
    tdat = deepcopy(s_tdat)
    vdat = deepcopy(s_vdat)
    indexfile = joinpath(filename, ".index")
    cf = pgc.readConfigFile(indexfile)
    if haskey(cf["."]["devices"], "Laser-Blue-raw")
        wavefunction =
            cf["."]["devices"]["Laser-Blue-raw"]["channels"]["pCell"]["waveGeneratorWidget"]["function"]
	    data_info["Laser.wavefunction"] = deepcopy(wavefunction)
    elseif haskey(cf["."]["devices"], "MultiClamp1")
	    wavefunction =
	        cf["."]["devices"]["MultiClamp1"]["waveGeneratorWidget"]["function"]
	    data_info["MultiClamp1.wavefunction"] = deepcopy(wavefunction)
		data_info["MultiClamp1.pulse_pars"] = [
			cf["."]["devices"]["MultiClamp1"]["waveGeneratorWidget"]["stimuli"]["Pulse"]["start"]["value"],
			cf["."]["devices"]["MultiClamp1"]["waveGeneratorWidget"]["stimuli"]["Pulse"]["length"]["value"],
			]
	
	else
        wavefunction = nothing
    end
    finalize(s_tdat)
    finalize(s_idat)
    finalize(s_vdat)
    finalize(wavefunction)
    @everywhere GC.gc()
    println("    Finished reading protocol data.")
    return tdat, idat, vdat, data_info
end

function get_lims(mode)
    #=
    Set the top and bottom limits to some defaults according to the 
    acquistion mode
    =#
    if mode == "'VC'"
        # print("VC")
        top_lims = (-0.4e-9, 0.4e-9)
        bot_lims = (-120e3, 100e3)
    elseif mode == "'IC'"
        top_lims = (-120e-3, 40e-3)
        bot_lims = (-2e-9, 2e-9)
    else
        println(RED_FG, "Unknown Mode: ", mode, WHITE_FG)
    end
end

function get_indices(data_info)
    #= 
    get the array indices that correpond to v and i
    depending on the acquisition mode
    =#
    mode = String(data_info["clampstate"]["mode"])  # will be VC, IC or IC=0 (may depend on age of acquisition code)
    c0_units = String(data_info["c0"]["units"])
    c1_units = String(data_info["c1"]["units"])
    if mode == "'VC'"
        # println("VCMode")
        if c0_units == "'A'"
            topidx = 1
            botidx = 2
        elseif c0_units == "'V'"
            topidx = 2
            botidx = 1
        end
    elseif mode == "'IC'"
        # println("ICMode")
        if c0_units == "'V'"
            topidx = 2
            botidx = 1
        elseif c0_units == "'A'"
            topidx = 1
            botidx = 2
        end
    else
        println(RED_FG, "mode is not known: ", mode, WHITE_FG)
    end
    return topidx, botidx
end

# function get_stim_times(data_info)
#     wv = data_info["Laser.wavefunction"]
#     u = split(wv, "\n")
#     stim_lats = Vector{Float64}()
#     re_float = r"[+-]?\d+\.?\d*"
#     for i = 1:size(u)[1]
#         s = match(re_float, u[i])
#         append!(stim_lats, parse(Float64, s.match) * 1e3)
#     end
#     return stim_lats
# end

function get_stim_times(data_info; device::AbstractString = "Laser")
    #=
    Retrieve stimulus times from a device's wavefunction parameters
    The default will get the information from a Laser device
    =#
    query = device * ".wavefunction"
    wv = data_info[query]
    u = split(wv, "\n")
    stim_lats = Vector{Float64}()
    re_float = r"[+-]?\d+\.?\d*"
    for i = 1:size(u)[1]
        s = match(re_float, u[i])
        append!(stim_lats, parse(Float64, s.match) * 1e3)
    end
    return stim_lats
end

function read_one_sweep(filename::AbstractString, sweep_dir, device)
    #=
    Get one waveform sweep from the protocol for
    the given sweep dir (protocol sequence) and the 
    specified device
    Returns the time array, the sweep data, and some info
    =#

    full_filename = joinpath(filename, sweep_dir, device)
    okfile = isfile(full_filename)
    if !okfile
        println(RED_FG, "File not found:")
        println("    ", full_filename, WHITE_FG)
        return false, false, false
    end
    tf = HDF5.ishdf5(full_filename)
    if !tf
        return false, false, false
    end
    fid = h5open(full_filename, "r")

    c0 = h5readattr(full_filename, "info/0/cols/0")
    c1 = h5readattr(full_filename, "info/0/cols/1")
    c2 = h5readattr(full_filename, "info/0/cols/1")
    ClampState = h5readattr(full_filename, "info/2/ClampState")
    DAQPrimary = h5readattr(full_filename, "info/2/DAQ/primary")
    DAQSecondary = h5readattr(full_filename, "info/2/DAQ/secondary")
    DAQCommand = h5readattr(full_filename, "info/2/DAQ/command")

    data_info = Dict(
        "clampstate" => deepcopy(ClampState),
        "c0" => c0,
        "c1" => c1,
        "c2" => c2,
        "DAQ.Primary" => deepcopy(DAQPrimary),
        "DAQ.Secondary" => deepcopy(DAQSecondary),
        "DAQ.Command" => deepcopy(DAQCommand),
        "Laser.wavefunction" => "",
    )

    time = deepcopy(fid["info"]["1"]["values"][:])
    data = deepcopy(fid["data"][:, :])
    close(fid)

    return time, data, data_info
end

function test_configread()
    filename = "/Users/pbmanis/Desktop/2018.09.27_000/ImageSequence_000/.index"
    data = pgc.readConfigFile(filename)
    println(data)
    fn = "/Users/pbmanis/Desktop/Python/mrk-nf107/data_for_testing/CCIV/.index"
    data = pgc.readConfigFile(fn)
    println(data)
end

function test_reader()
    # local test file name
    #filename = "/Volumes/Pegasus/ManisLab_Data3/Kasten_Michael/Pyramidal/2018.02.12_000/slice_001/cell_000/CCIV_1nA_max_000"
    filename = "/Users/pbmanis/Desktop/Python/mrk-nf107/data_for_testing/CCIV_1nA_max_1s_pulse_000"
    read_hdf5(filename)
end


end
