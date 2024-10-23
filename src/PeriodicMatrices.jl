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
#function  pseigsm end
#function  psceigsm end
#function  psm2hr end
#function  hr2psm end

"""
    psceigfr(A::FourierFunctionMatrix[, N]; P, atol) -> ce

Compute the characteristic exponents of a continuous-time periodic matrix in _Fourier Function Matrix_ representation. 

For a given square continuous-time periodic function matrix `A(t)` of period `T` 
in a  _Fourier Function Matrix_ representation, 
the characteristic exponents `ce` are computed as the eigenvalues of the state operator `A(t)-D*I` lying in the 
fundamental strip `-ω/2 <  Im(λ) ≤ ω/2`, where `ω = 2π/T`. A finite dimensional truncated matrix of order `n*(2*N*P+1)` 
is built to approximate `A(t)-D*I`, where `n` is the order of `A(t)`,  `N` is the number of selected harmonic components
in the Fourier representation and `P` is the period multiplication number (default: `P = 1`).
The default value used for `N` is `N = max(10,p-1)`, where `p` the number of harmonics terms of `A(t)` (see [`FourierFunctionMatrix`](@ref)). 

The keyword argument `atol` (default: `atol = 1.e-10`) is a tolerance on the magnitude of the trailing components of the 
associated eigenvectors used to validate their asymptotic (exponential) decay. Only eigenvalues satisfying this check are returned in `ce`. 
"""
function  psceigfr end
"""
     ffm2psm(Af::FourierFunctionMatrix, nrange atol = 0, rtol = 10*n*ϵ,) -> A::Matrix{Num}

Convert a range of harmonic components `nrange` of the Fourier function matrix `Af` to a symbolic matrix `A`. 
The default range is `nrange = 0:n`, where `n` is the order of the maximum harmonics. 
The tolerance used to assess nonzero coefficients is `tol = max(atol, rtol*maxnorm)`, where 
`maxnorm` is the maximum absolute value of the coefficients `ac_i` and `as_i` in `Af(t)`. The default values of tolerances
are `atol = 0` and `rtol = 10*n*ϵ`, where `ϵ` is the working machine precision.
"""
function  ffm2psm end
"""
     ffm2hr(Afun::FourierFunctionMatrix; atol = 0, rtol = √ϵ, squeeze = true) -> Ahr::HarmonicArray

Compute the harmonic (Fourier) representation of a Fourier periodic matrix object. 

The Fourier function matrix object `Afun` of period `T` is built using
the Fourier series representation of a periodic matrix `Afun(t)` of subperiod `T′ = T/k`, 
where each entry of `Afun(t)` has the form

             p
      a_0 +  ∑ ( ac_i*cos(i*t*2*π/T′)+as_i*sin(i*2*π*t/T′) ) ,
            i=1 

where `k ≥ 1` is the number of subperiods (default: `k = 1`).   

The harmonic array object `Ahr` of period `T` is built using
the harmonic representation of a periodic matrix `Ahr(t)` of subperiod `T′′ = T/k′` in the form

                     p′
     Ahr(t) = A_0 +  ∑ ( Ac_i*cos(i*t*2*π/T′′)+As_i*sin(i*2*π*t/T′′) ) ,
                    i=1 

where `p′` is the maximum index `j`, such that `Ac_j` and/or `As_j` are nonzero.
The tolerance used to assess nonzero elements is `tol = max(atol, rtol*maxnorm)`, where 
`maxnorm` is the maximum absolute value of the coefficients `ac_i` and `as_i` in `Afun(t)`. The default values of tolerances
are `atol = 0` and `rtol = √ϵ`, where `ϵ` is the working machine precision.
The resulting harmonic approximation `Ahr(t)` is returned in the harmonic array object `Ahr` 
(see [`HarmonicArray`](@ref)). 
"""
function  ffm2hr end



include("types/PeriodicMatrix.jl")
include("types/PeriodicMatrix_sym.jl")
#include("types/PeriodicMatrix_Fourier.jl")
#include("SymFourierConv.jl")
include("pmops.jl")
include("pmops_sym.jl")
#include("pmops_Fourier.jl")
include("psfutils.jl")
include("pmutils.jl")
include("pmutils_sym.jl")
#include("pmutils_Fourier.jl")
include("pmutilities.jl")

end
