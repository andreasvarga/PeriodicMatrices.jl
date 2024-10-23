module FourierApproxExt

using PeriodicMatrices
using ApproxFun
using LinearAlgebra
using Symbolics
using OrdinaryDiffEq
using IRKGaussLegendre
using Primes
using QuadGK
using Optim


import Base: +, -, *, /, \, (==), (!=), ^, isapprox, iszero, isequal, convert, promote_op, size, length, ndims, reverse, 
             hcat, vcat, hvcat, inv, show, lastindex, require_one_based_indexing, print, show, one, zero, eltype

export FourierFunctionMatrix, psceigfr, ffm2hr

include("types/PeriodicMatrix_Fourier.jl")
include("pmops_Fourier.jl")
include("pmutils_Fourier.jl")
include("SymFourierConv.jl")

end
