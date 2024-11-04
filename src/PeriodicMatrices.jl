module PeriodicMatrices

# using ApproxFun
using LinearAlgebra
using FFTW
using Interpolations
using Optim
using OrdinaryDiffEq
using IRKGaussLegendre
using Primes
using QuadGK
using PeriodicSchurDecompositions
using SLICOT_jll
using Symbolics

include("SLICOTtools.jl")
using .SLICOTtools: mb03vd!, mb03vy!, mb03bd!, mb03wd!, mb03vw!, mb03kd!


import LinearAlgebra: BlasInt, BlasFloat, BlasReal, BlasComplex, copy_oftype, transpose, adjoint, opnorm, normalize, rdiv!, issymmetric, norm, tr
import Base: +, -, *, /, \, (==), (!=), ^, isapprox, iszero, isequal, convert, promote_op, size, length, ndims, reverse, 
             hcat, vcat, hvcat, inv, show, lastindex, require_one_based_indexing, print, show, one, zero, eltype

export PeriodicMatrix, pschur, pschur!, pschur1, pschur2, pgschur, pgschur!, phess, phess!, phess1, psreduc_reg, psreduc_fast, check_psim, mshift,  
       tvmeval, tpmeval, hreval, tvstm, psordschur!, psordschur1!, pgordschur!
export ts2hr, ts2pfm, ts2fm, tsw2pfm, ts2ffm, pfm2hr, pm2pa, ffm2hr, pmaverage, hrtrunc, hrchop
export monodromy, pseig, psceig, psceighr, psceigfr, peigvals
export PeriodicArray, PeriodicMatrix, SwitchingPeriodicArray, SwitchingPeriodicMatrix
export PeriodicTimeSeriesMatrix, PeriodicSwitchingMatrix, HarmonicArray, PeriodicFunctionMatrix
export isconstant, iscontinuous, isdiscrete, set_period, promote_period, promote_period2, getpm
export mb03vd!, mb03vy!, mb03bd!, mb03wd!, mb03kd! 
export ps2fls, hr2bt, hr2btupd, phasemat, ps2frls, DiagDerOp, ps2ls, ps2spls
export pmshift, trace, pmzeros, pmcopy
export pmderiv, pmrand, horzcat, vertcat, pmsymadd!
export bldiag, blockdiag, blockut
export pmmulsym, pmtrmulsym, pmmultrsym, pmmuladdsym, pmmultraddsym, pmmuladdtrsym, pmata, pmaat
export AbstractPeriodicArray

export PeriodicSymbolicMatrix, pseigsm, psceigsm, psm2hr, hr2psm
export FourierFunctionMatrix, psceigfr, ffm2hr, ffm2psm

abstract type AbstractPeriodicArray{Domain,T} end

include("types/PeriodicMatrix.jl")
include("types/PeriodicMatrix_sym.jl")
include("pmops.jl")
include("pmops_sym.jl")
include("psfutils.jl")
include("pmutils.jl")
include("pmutils_sym.jl")
include("pmutilities.jl")
include("docstrings.jl")
function DiagDerOp end

end
