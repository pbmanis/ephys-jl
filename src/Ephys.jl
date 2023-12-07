using Ephys

push!(ARGS, "arg")
Acq4Reader.julia_main()

# It is ok to use stdlibs that are not in the project dependencies
using Test
@test 1==1