module Configfile
# using PyCall
# using JSON
using Base.Meta

export read_configFile, filetest, test
# configfile.jl : based on:
# 	configfile.py - Human-readable text configuration file library
# 	Copyright 2010  Luke Campagnola
# 	Distributed under MIT/X11 license. See license.txt for more infomation.
#
# Used for reading a python-like configuration
# file format, as created and used in acq4. Data structures may be nested and contain any data type as long
# as it can be converted to/from a string using repr and eval.
#
# The configfile is just a text file.
# This is a pretty straight translation of the original code.

## Very simple unit support:
##  - creates variable names like 'mV' and 'kHz'
##  - the value assigned to the variable corresponds to the scale prefix
##    (mV = 0.001)
##  - the actual units are purely cosmetic for making code clearer:
##  
##    x = 20*pA    is identical to    x = 20*1e-12

## No unicode variable names (μ,Ω) allowed until python 3

SI_PREFIXES = "yzafpnum kMGTPEZY"
UNITS = split("m,s,g,W,J,V,A,F,T,Hz,Ohm,S,N,C,px,b,B", ",")
allUnits = Vector{String}()
cm = 0.01

function addUnit(p, n)
    v = 1000^n
    for u in UNITS
        # g[p+u] = v
        allUnits[p + u] = v
    end
end

function make_SI_array(p)
    for p in SI_PREFIXES
        if p === " "
            p = ""
            n = 0
        elseif p === "u"
            n = -2
        else
            n = SI_PREFIXES.index(p) - 8
        end
        addUnit(p, n)
    end
end

function evalUnits(unitStr)
    """
    Evaluate a unit string into ([numerators,...], [denominators,...])
    Examples:
        N m/s^2   =>  ([N, m], [s, s])
        A*s / V   =>  ([A, s], [V,])
    """
    return pass
end

function formatUnits(units)
    """
    Format a unit specification ([numerators,...], [denominators,...])
    into a string (this is the inverse of evalUnits)
    """
    return pass
end

function simplify(units)
    """
    Cancel units that appear in both numerator and denominator, then attempt to replace 
    groups of units with single units where possible (ie, J/s => W)
    """
    return pass
end

function ColorMap(c)
    c = replace(c, "ColorMap(array([" => "([")
    c = replace(c, ", dtype=uint8)" => "")
    (p, c) = split(c, ", array("; limit=2)
    c = replace(c, "array(" => "")
    c = replace(c, ")" => "")
    return p, c
end

function Point(p)
    return list(p)
end

function measureIndent(s::AbstractString)
    n = 1
    while n <= length(s) && s[n] == ' '
        n += 1  # count the number of spaces from start of line
    end
    # println("length s: ", length(s), "  n = ", n, "  s=<", s, ">")
    if length(n) < length(s) && s[n] == '#' # comment line only
        return 0
    end
    if length(s) == n - 1
        return 0
    end # blank line
    return n - 1
end

function parseConfigString(lines::AbstractString; start::Int64=1, verbose::Bool=false)
    data = Dict()  # store parsed data at this level in a dictionary

    splitlines = [li for li in split(lines, "\n"; keepempty=false)] # break by line breaks
    splitlines = [
        sl for (i, sl) in enumerate(splitlines) if
        occursin(r"\S", sl) & !(match(r"\s*(?:#|$)", sl) === nothing)
    ]  ## remove empty lines and comments

    indent = measureIndent(splitlines[start])  # get the current indent level
    ln = start - 1  # because we will increment this in the while true ... 

    while true
        ln += 1
        if ln >= length(splitlines)
            break
        end
        thisline = splitlines[ln]

        ## Skip blank lines or lines starting with #
        if !occursin(r"\S", thisline)
            continue
        end

        ## Measure line indentation, make sure it is correct for this level
        lineInd = measureIndent(thisline)

        if lineInd < indent
            ln -= 1 # back up one level
            break  # will cause return
        end
        if lineInd > indent
            println(
                "    Indentation is incorrect. Expected ",
                indent,
                "  got: ",
                lineInd,
                "  on line: ",
                ln,
            )
            println("        Offending string: ", thisline)
            error("Bad indent")
        end
        if length(thisline) == 0 | isnothing(thisline)
            continue
        end
        # check for comment and remove it
        if !isnothing(match(r"[#]+", thisline))
            thisline = split(thisline, "#"; limit=2)[1]  # strip off the comment at the end of the line
            if !isnothing(match(r"[\s]+", thisline))
                thisline = nothing
            end
        end

        if !isnothing(thisline)
            # check for key
            if isnothing(match(r"[:]+", thisline))
                println("    Missing colon on line: ", ln)
                println("    Offending string: ", thisline)
                error("Missing colon")
            end
            values = split(thisline, ':'; limit=2)
            values = strip.(values)
            key = values[1]
            if length(key) == 0
                error("Missing name preceding colon")
            end
            # get the right-hand side of the key and re-interpret it
            value = values[2]
            value = replace(value, "\'" => "\"")  # replace ' with "
            value = replace(value, "False" => "false")
            value = replace(value, "True" => "true")
            value = replace(value, "u\"None\"" => "None")
            if (length(value) > 0)  # must have a value
                # remove "L" ("long") from the values
                if (value[1] in ['[', '(']) & (value[end] in [']', ")"])
                    value = replace(value, "L" => "")
                end
                # strip out the "array" as as Meta.parse cannot evaluate it
                # also convert all "None" to "nothing" for julia
                if (value[1] in ['{', '(', '[']) & (value[end] in [']', ')', '}'])
                    value = replace(value, "array([" => "([")
                    value = replace(value, "None" => "nothing")
                    # println("Array value to parse: ", value)
                    #val = py"eval_value"(value)
                    val = eval(Meta.parse(value))
                elseif startswith(value, "ColorMap")
                    p, value = ColorMap(value)
                    value = join([p, value], ", ")
                    val = eval(Meta.parse(value))
                elseif startswith(value, "#")
                    val = nothing
                else
                    # val = py"eval_value"(value)
                    val = string(value) # eval(Meta.parse(value))
                end

            else  # no argument/value
                # if key == "."
                #     val = nothing
                #     println(" setting root key ")
                if (ln + 1 > length(lines)) | (measureIndent(splitlines[ln + 1]) <= indent)
                    # println( "    blank dict for key: ", key)
                    val = nothing
                else
                    # println(" >>>> Going deeper..", ln)
                    ln, val = parseConfigString(lines; start=ln + 1)
                end
            end
            # println("setting key: ", key, " to value: ", val)
            data[key] = val
            # println("\n\nKEYS:\n", keys(data))
        end
    end
    return ln, data
end

function read_configfile(fname)
    try
        open(fname, "r") do fd
            s = read(fd, String)
            s = replace(s, "\r\n" => s"\n")
            s = replace(s, "\r" => s"\n")
            # println("Calling parseString")
            ln, data = parseConfigString(s)
            # println("\n\n")
            # println(keys(data))

            return data
        end
    catch
        println("Error while reading config file: ", fname)
        error()
    end
end

# Not Implemented:
# def writeConfigFile(data, fname):
#     s = genString(data)
#     fd = open(fname, 'w')
#     fd.write(s)
#     fd.close()
#
# def appendConfigFile(data, fname):
#     s = genString(data)
#     fd = open(fname, 'a')
#     fd.write(s)
#     fd.close()
#
#
# def genString(data, indent=''):
#     s = ''
#     for k in data:
#         sk = str(k)
#         if len(sk) == 0:
#             print(data)
#             raise Exception('blank dict keys not allowed (see data above)')
#         if sk[0] == ' ' or ':' in sk:
#             print(data)
#             raise Exception('dict keys must not contain ":" or start with spaces [offending key is "%s"]' % sk)
#         if isinstance(data[k], dict):
#             s += indent + sk + ':\n'
#             s += genString(data[k], indent + '    ')
#         else:
#             s += indent + sk + ': ' + repr(data[k]) + '\n'
#     return s
#

function test()
    fn = "/Volumes/Pegasus_002/ManisLab_Data3/Kasten_Michael/NF107Ai32_Het/2018.02.12_000/slice_001/cell_000/CCIV_4nA_max_002/.index"
    return data = read_configfile(fn)
    # print(json(data, 4))
    # fn = "/Users/pbmanis/Desktop/Python/mrk-nf107/data_for_testing/CCIV/.index"
    # data = read_configfile(fn)
    # Println(data)
end

function filetest()
    fn, fnio = mktemp()

    cf = """
    key2:              #comment
                    #comment
    key21: 'value' # comment
                    #comment
    key22: [1,2,3]
    key23: 234  #comment
    key99: 'value'
    """
    write(fnio, cf)
    close(fnio)
    println("=== Test:===")
    num = 0
    for line in eachsplit(cf, "\n")
        num += 1
        println(num, ": ", line)
    end
    println(cf)
    println("============")

    return data = read_configfile(fn)
end

end
