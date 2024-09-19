module Runtests

using Test
using PeriodicMatrices

@testset "Test PeriodicMatrices" begin
# test constructors, basic tools
include("test_pschur.jl")
include("test_pmutils.jl")
include("test_pmops.jl")
end

end
