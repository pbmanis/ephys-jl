using Acq4Reader
using Test

@testset "Acq4Reader.jl" begin
    # Write your tests here.
    # test 1: did module load ok?
    # err = Acq4Reader.test_configread()
    err = Acq4Reader.test_reader()
#@test err == false

end
