module Acq4Reader
using Statistics
using HDF5
using Plots
using Printf

using InteractiveUtils

export read_hdf5, get_lims


pyplot()
fname = ""


function greeter()
    println("Acq4Reader_Module")
end

greeter()

device = "MultiClamp1.ma"

# filename = "/Volumes/Pegasus_002/Additional_Kasten/2020.02.24_000/slice_001/cell_000/fast_IO_CC_001/"
# filename = "/Volumes/Pegasus_002/Additional_Kasten/2020.02.24_000/slice_000/cell_000/fast_IO_001/"
# sweep_dir = "000_000"
filename = "/Volumes/Pegasus_002/ManisLab_Data3/Kasten_Michael/ChATIREScre_ai32/2021.02.17_000/slice_000/cell_001/CCIV_long_001/"


function get_subdirs(base_dir)
    # read the subdirectories of the current base_dir
    # only return those that are protocol sweeps
    subdirs = readdir(base_dir)
    deleteat!(subdirs, subdirs .== ".index")
    return subdirs
end
#
# sweeps = get_subdirs(filename)

function read_hdf5(filename)
    #=
    Read the protocol data from the specified metaarray file,
    treating it as an HDF5 file (which it is)
        
    Return the data, and some "standard" y limits.
    =#
    sweeps = get_subdirs(filename)
    idat = Array{Float64}(undef)
    vdat = Array{Float64}(undef)
    tdat = Array{Float64}(undef)
    first = true
    mode = ""
    clampstate = ""
    data_info = ""
    @time for s in sweeps
        time, data, data_info = read_one_sweep(filename, s, device)
        if time == false
            continue
        end
        sweep_mode = String(data_info["clampstate"]["mode"])  # will be VC, IC or IC=0 (may depend on age of acquisition code)
        indx1, indx2 = get_indices(data_info)

        if first
            first = false
            tdat = time
            idat = data[:, indx1]
            vdat = data[:, indx2]
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
            tdat = hcat(tdat, time)
            idat = hcat(idat, data[:, indx1])
            vdat = hcat(vdat, data[:, indx2])
        end
    end
    if mode == ""
        throw(ErrorException("No acquisition mode found"))
    end
    println("Protocol reading finished.")
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
        println("Unknown Mode HELP")
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
        println("mode is not known: ", mode)
    end
    return topidx, botidx
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
        println("File not found!!!!!")
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

    data_info =
        Dict("clampstate" => deepcopy(ClampState), "c0" => c0, "c1" => c1, "c2" => c2)

    time = fid["info"]["1"]["values"][:]
    data = deepcopy(fid["data"][:, :])
    close(fid)
    return time, data, data_info
end




end
