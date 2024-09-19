```@meta
CurrentModule = PeriodicMatrices
DocTestSetup = quote
    using PeriodicMatrices
end
```

# PeriodicMatrices.jl

[![DocBuild](https://github.com/andreasvarga/PeriodicMatrices.jl/workflows/CI/badge.svg)](https://github.com/andreasvarga/PeriodicMatrices.jl/actions)
[![Code on Github.](https://img.shields.io/badge/code%20on-github-blue.svg)](https://github.com/andreasvarga/PeriodicMatrices.jl)

`PeriodicMatrices.jl` provides the basic tools to handle periodic time-varying matrices. 
The time dependence can be either continuous or discrete. 

Periodic matrices appear in many control applications involving periodic phenomena as for example, satellite attitude control, helicopter forward flight control, 
orbital stabilization of underactuated systems, control of multi-rate sampled-data systems, etc. These problems can be addressed using tools available in 
the [PeriodicSystems](https://github.com/andreasvarga/PeriodicSystems.jl) package, which is based on the periodic matrix type defined in this package.  

A continuous-time periodic matrix can be specified in the following forms:

- periodic matrix function
- harmonic matrix series
- periodic matrix time series with uniform time grid 
- periodic matrix time series with non-uniform time grid (also known as periodic switching matrix)
- periodic symbolic matrix
- Fourier matrix series approximation   

A discrete-time periodic matrix can be specified in the following forms:

- periodic matrix time series with time-varying dimensions with uniform time grid
- periodic matrix time series with time-varying dimensions with non-uniform time grid
- periodic matrix time series with constant dimensions with uniform time grid
- periodic matrix time series with constant dimensions with non-uniform time grid

For a periodic matrix `A(t)` of period `T` it is not assumed that `T` is the minimum value
which satisfies the periodicity condition `A(t) = A(t+T)` for all values of `t`. To describe 
matrices having multiple periods, a subperiod `Tsub := T/n` can be defined, such that `A(t) = A(t+Tsub)`,
for all `t`. This allows a substantial memory saving for some classes of periodic representations. 

Several operations on periodic matrices are implemented, such as, inversion, transposing, norms, derivative/shifting, trace.
All operations with two periodic matrices such as addition/substraction, multiplication, horizontal/vertical concatenation, block-diagonal appending,
allow different, but commensurate, periods/subperiods.  

Functions are provided to compute the characteristic multipliers and characteristic exponents of periodic matrices, using methods based on the periodic Schur decomposition of matrix products 
or structure exploitung fast algorithms. 
These functions are instrumental to apply [Floquet theory](https://en.wikipedia.org/wiki/Floquet_theory) to study the properties of solutions of 
various classes of differential equations (Mathieu, Hill, Meissner) and the stability of linear periodic systems (see [PeriodicSystems](https://github.com/andreasvarga/PeriodicSystems.jl) package). 


The following categories of functions are currently implemented:

**Constructors for periodic matrices**

* **[`PeriodicMatrix`](@ref)**   Discrete-time periodic matrix representation.
* **[`PeriodicArray`](@ref)**    Discrete-time periodic array representation.
* **[`SwitchingPeriodicMatrix`](@ref)** Discrete-time switching periodic matrix representation.
* **[`SwitchingPeriodicArray`](@ref)** Discrete-time switching periodic array representation.
* **[`PeriodicFunctionMatrix`](@ref)**  Continuous-time periodic function matrix representation.
* **[`PeriodicSymbolicMatrix`](@ref)**   Continuous-time periodic symbolic matrix representation.
* **[`PeriodicTimeSeriesMatrix`](@ref)**   Continuous-time periodic time series matrix representation.
* **[`HarmonicArray`](@ref)**   Continuous-time harmonic array representation.
* **[`FourierFunctionMatrix`](@ref)**   Continuous-time Fourier functin matrix representation.
* **[`PeriodicSwitchingMatrix`](@ref)** Continuous-time switching periodic matrix representation.

**Periodic matrix conversions**

* **[`ts2hr`](@ref)**   Conversion of  a periodic time series matrix to a harmonic array approximation.
* **[`pfm2hr`](@ref)**  Conversion of  a periodic function matrix to a harmonic array representation. 
* **[`ts2pfm`](@ref)**  Conversion of  an interpolated periodic time series matrix to a periodic function matrix.
* **[`hr2psm`](@ref)**  Conversion of  a harmonic array representation to a periodic symbolic matrix.
* **[`psm2hr`](@ref)**  Conversion of  a periodic symbolic matrix into a harmonic array representation.
* **[`pm2pa`](@ref)**   Conversion of  a discrete-time periodic matrix object to a periodic array object.
* **[`ffm2hr`](@ref)**  Conversion of  a Fourier function matrix to a harmonic array representation. 
* **[`hr2bt`](@ref)**   Building a block Toeplitz matrix approximation of a harmonic (Fourier) array representation. 
* **[`hr2btupd`](@ref)**  Building an updated block Toeplitz matrix approximation of a harmonic (Fourier) array representation. 

**Periodic matrix utilities**

* **[`pseig`](@ref)**   Characteristic multipliers of a periodic matrix.
* **[`psceig`](@ref)**   Characteristic exponents of a periodic matrix.
* **[`psceighr`](@ref)**   Characteristic exponents of a periodic matrix in Harmonic Array representation.
* **[`psceigfr`](@ref)**   Characteristic exponents of a periodic matrix in Fourier Function Matrix representation.
* **[`monodromy`](@ref)**  Monodromy matrix of a linear periodic time-varying system of ODE.
* **[`tvstm`](@ref)**  State transition matrix of a linear time-varying system of ODE.
* **[`psreduc_reg`](@ref)**  Fast reduction of a lifted regular pencil corresponding to a product of matrices. 
* **[`tvmeval`](@ref)**  Time response evaluation of a continuous-time periodic matrix. 
* **[`hreval`](@ref)**  Evaluation of a harmonic array for a numerical or symbolic time value. 
* **[`hrchop`](@ref)**  Removal of the negligible trailing terms of a harmonic representation. 
* **[`hrtrunc`](@ref)**  Truncation of a harmonic representation.  
* **[`pmaverage`](@ref)**  Evaluation of the time averaged matrix of a continuous-time periodic matrix. 

**Periodic Schur decompositions**

* **[`phess`](@ref)**  Periodic Hessenberg decomposition of a product of matrices.
* **[`pschur`](@ref)**  Periodic Schur decompositions of products or quotient products of matrices. 
* **[`psordschur!`](@ref)**  Reordering of periodic Schur decompositions of products or quotient products of matrices.
* **[`psordschur1!`](@ref)**  Reordering of periodic Schur decompositions of products or quotient products of square matrices.
* **[`pgschur`](@ref)**  Generalized real periodic Schur decomposition of a formal product of matrices.
* **[`pgschur`](@ref)**  Generalized real periodic Schur decomposition of a formal product of matrices.
* **[`pgschur!`](@ref)**  Generalized real periodic Schur decompositions of formal products of matrices (in place computation).
* **[`pgschur`](@ref)**  Generalized real periodic Schur decompositions of formal products of matrices.
* **[`pgordschur!`](@ref)**  Reordering of generalized real periodic Schur decompositions a formal products of matrices.

## [Release Notes](https://github.com/andreasvarga/PeriodicMatrices.jl/blob/master/ReleaseNotes.md)

## Main developer

[Andreas Varga](https://sites.google.com/view/andreasvarga/home)

License: MIT (expat)

## References

[1] A. Varga. [A Periodic Systems Toolbox for Matlab](https://elib.dlr.de/12283/1/varga_ifac2005p1.pdf). Proc. of IFAC 2005 World Congress, Prague, Czech Republic, 2005.

[2] S. Bittanti and P. Colaneri. Periodic Systems - Filtering and Control, Springer Verlag, 2009.

[3] J. A. Richards. Analysis of Periodically Time-Varying Systems, Springer Verlag, 1983.
