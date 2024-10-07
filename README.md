# PeriodicMatrices.jl

<!-- [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.4568159.svg)](https://doi.org/10.5281/zenodo.4568159) -->
[![codecov.io](https://codecov.io/gh/andreasvarga/PeriodicMatrices.jl/coverage.svg?branch=main)](https://codecov.io/gh/andreasvarga/PeriodicMatrices.jl?branch=main)
[![Latest](https://img.shields.io/badge/docs-latest-blue.svg)](https://andreasvarga.github.io/PeriodicMatrices.jl/dev/)
[![The MIT License](https://img.shields.io/badge/license-MIT-brightgreen.svg?style=flat-square)](https://github.com/andreasvarga/PeriodicMatrices.jl/blob/main/LICENSE.md)
[![CI](https://github.com/andreasvarga/PeriodicMatrices/actions/workflows/CI.yml/badge.svg)](https://github.com/andreasvarga/PeriodicMatrices/actions/workflows/CI.yml)

## Handling of periodic time-varying matrices

## Compatibility

Julia 1.10 and higher.

<!-- ## How to install

````JULIA
pkg> add PeriodicMatrices
pkg> test PeriodicMatrices
```` -->

## About

`PeriodicMatrices.jl` provides the basic tools to handle periodic time-varying matrices. 
The time dependence can be either continuous or discrete. 

A continuous-time periodic matrix can be specified in the following forms:

- periodic matrix function
- harmonic matrix series
- periodic matrix time series with uniform time grid 
- periodic matrix time series with non-uniform time grid
- periodic symbolic matrix
- Fourier matrix series approximation   

A discrete-time periodic matrix can be specified in the following forms:

- periodic matrix time series with time-varying dimensions with uniform time grid
- periodic matrix time series with time-varying dimensions with non-uniform time grid
- periodic matrix time series with constant dimensions with uniform time grid
- periodic matrix time series with constant dimensions with non-uniform time grid

For a periodic matrix `A(t)` of period `T` it is _not_ assumed that `T` is the minimum value
which satisfies the periodicity condition `A(t) = A(t+T)` for all values of `t`. To describe 
matrices having multiple periods, a subperiod `Tsub := T/n` can be defined, such that `A(t) = A(t+Tsub)`,
for all `t`. This allows a substantial memory saving for some classes of periodic representations. 

The provided classes of periodic representation extend the classes used in the _Periodic Systems Toolbox for Matlab_ (see [1]).  

Several operations on periodic matrices are implemented, such as, inversion, transposing, norms, derivative/shifting, trace.
All operations with two periodic matrices such as addition/substraction, multiplication, horizontal/vertical concatenation, block-diagonal appending,
allow different, but commensurate, periods/subperiods.  

Functions are provided to compute the characteristic multipliers and characteristic exponents of periodic matrices, using methods based on the periodic Schur decomposition of matrix products 
or structure exploitung fast algorithms. 
These functions are instrumental to apply [Floquet theory](https://en.wikipedia.org/wiki/Floquet_theory) to study the properties of solutions of 
various classes of differential equations (Mathieu, Hill, Meissner) and the stability of linear periodic systems (see [PeriodicSystems](https://github.com/andreasvarga/PeriodicSystems.jl) package). 
 
## Example: Floquet-analysis of differential equations with periodic parameters

A frequently encountered periodic differential equation is of second order, expressible as

$$\ddot{x} + 2\zeta\dot{x} + (a - 2q\psi(t))x = 0 ,$$

where  $ψ(t) = ψ(t+T)$ is a periodic function of period $T$. The parameter $a$ represents a constant portion of the
coefficient of $x$ and $q$ accounts for the magnitude of the time variation. A positive damping coefficient $\zeta$ is frequently used in practical applications.  
The above equation for $\zeta = 0$ is generally known as the [_Hill-equation_](https://en.wikipedia.org/wiki/Hill_differential_equation) and the form in which
it is expressed is that most widely encountered in applications. When $\psi(t) = cos 2t$ , the equation becomes the [_Mathieu equation_](https://en.wikipedia.org/wiki/Mathieu_function#Mathieu_equation).
If $\phi(t)$ is a rectangular function, the corresponding particular form is known
as the [_Meissner equation_](https://en.wikipedia.org/wiki/Meissner_equation). If $\zeta > 0$, we speak of the _lossy_ variants of these equations.

The above equation can be equivalently expressed as a first order system of differential equations, by defining

```math
y(t) = \left[  \begin{array}{c} x(t)\\ \dot{x}(t) \end{array} \right] 
```    

to recast the second order differential equation into the form

$$ \dot{y}(t) = A(t)y(t)$$

where $A(t)$ is a $2\times 2$ matrix of constant or periodically varying coefficients of the form

```math
   A(t) = \left[ \begin{array}{cc} 
          0 & 1\\
          -a+2q\psi(t) & -2\zeta 
          \end{array} \right]
```         

The state transition matrix $\Phi(t,0)$ over the time interval $(0,t)$ satisfies the differential equation 

```math
\dot{\Phi}(t,0) = A(t)\Phi(t,0),  \,\, \Phi(0,0) = I 
``` 

and the _monodromy matrix_ $\Psi := \Phi(T,0)$, i.e., the state transition matrix over one full period.

The Floquet-theory based stability analysis of the above equations addresses the determination of _characteristic multipliers_ $\lambda_i$ as the eigenvalues of the monodromy matrix
or alternatively the _characteristic exponents_ $\mu_i$ related to the characteristic multipliers as

$$ \lambda_i = exp(\mu_iT) .$$ 

The solution $x(t)$ is _stable_ if it remains bounded as time goes to infinity.  For stability it is sufficient that $Re(\mu_i) < 0 \, \forall i$, which is the same as $|\lambda_i| < 1 \, \forall i$. 
Such a solution will also be _stable_ if in addition one $Re(\mu_i) = 0$ or one $|\lambda_i| = 1$. 

In what follows we illustrate how to perform the stability analysis for the three types of equations based on the characteristic multipliers/exponents.  

### The lossless Meissner equation with a rectangular waveform coefficient

This is the Example of Fig 3.1 in [3]. 
Assume the period is $T = \pi$ and let $\tau = T/3$ be the switching time. We consider the periodic function $\psi(t) = 1$ if $t \in [0,\tau)$ and $\psi(t) = -1$ if $t \in [\tau,\pi)$. 
We can describe the periodic matrix $A(t)$ as a **PeriodicSwitchingMatrix**  with two components corresponding to the two constant values of $\psi(t)$ and switching times at $t = 0$ and $t = \tau$. 
The following code can be used for stability analysis purposes:

````JULIA
using PeriodicMatrices

# set parameters
a = 1; q = .1; T = pi
ts = [0; τ; T] 

ψ(t,ts) = isodd(findfirst(mod(t,T) .< ts)) ? -1 : 1

# setup of A(t)
A = PeriodicSwitchingMatrix([[0. 1.; -a+2*q*ψ(t,ts) 0] for t in ts[1:end-1]], ts[1:end-1], T)
ce = psceig(A)  

# stability test
all(real(ce) .< 0)
````

The computed characteristic exponents are:

````JULIA
julia> ce 
2-element Vector{ComplexF64}:
 0.041379744661220644 + 1.0im
  -0.0413797446612206 + 1.0im
````
and therefore the solutions are unstable. The computations can be easily extended to several switching points in $\tau$ as well.

For a lossy Meissner equation with $\zeta = 0.2$, the computations reveal stability:

````JULIA
# setup of A(t)
ζ = 0.2
A = PeriodicSwitchingMatrix([[0. 1.; -a+2*q*ψ(t,ts) -2ζ] for t in ts[1:end-1]], ts[1:end-1], T)
ce = psceig(A)  

# stability test
all(real(ce) .< 0)
````
and 
````JULIA
julia> ce
2-element Vector{ComplexF64}:
 -0.14798307933724503 + 1.0im
   -0.252016920662755 + 1.0im
````

### The Hill equation with a sawtooth waveform coefficient

This is the Example of Fig 3.3 in [3]. Assume again the period is $T = \pi$ and let $\psi(t)$ be the periodic negative slope sawtooth function $\psi(t) = -2*t/T+1$. 
We can describe the periodic matrix $A(t)$ as a **PeriodicFunctionMatrix**. 
The following code can be used for stability analysis purposes:

````JULIA
using PeriodicMatrices

# set parameters
a = 1; q = .1; T = pi
ψ(t) = -2*t/pi+1

# setup of A(t)
A = PeriodicFunctionMatrix(t->[0. 1.; -a+2*q*ψ(t) 0], T)
ce = psceig(A)  

# stability test
all(real(ce) .< 0)
````
The computed characteristic exponents are:

````JULIA
julia> ce
2-element Vector{ComplexF64}:
 0.031815620244098085 + 1.0im
  -0.0318155972239221 + 1.0im
````
and therefore the solutions are unstable. 

Stability can be illustrated with the lossy Hill equation $\zeta = 0.2$:  

````JULIA
# setup of A(t)
ζ = 0.2
A = PeriodicFunctionMatrix(t->[0. 1.; -a+2*q*ψ(t) -2ζ], T)

ce = psceig(A)  

# stability test
all(real(ce) .< 0)
````
and 
````JULIA
ulia> ce = psceig(A)
2-element Vector{ComplexF64}:
 -0.17487078757796817 + 1.0im
 -0.22513738474355782 + 1.0im
````

Alternatively, the same computations can be performed with $A$ defined as a **PeriodicSymbolicMatrix**: 

````JULIA
using PeriodicMatrices
using Symbolics

# set parameters
a = 1; q = .1; T = pi
@variables t
ψ1 = -2*t/pi+1

# setup of A(t)
A = PeriodicSymbolicMatrix([0. 1.; -a+2*q*ψ1 0], T)
ce = psceig(A)  

# stability test
all(real(ce) .< 0)
````

### Using multiple shooting to enhance accuracy

Consider the following periodic matrix of period $T = 2pi$


```math
   A(t) = \left[ \begin{array}{cc} 
          0 & 1\\
          -10*cos(t) & -24-10*sin(t) 
          \end{array} \right]
```         
which has characterisitic exponents equal to $0$ and $-24$.

Using the standard settings to compute the characteristic exponents,  the resulting characteristic exponents have no one exact digit 
even by imposing a high relative tolerance for solving the underlying differential equations:

````JULIA
julia> using Periodicmatrices

julia> A = PeriodicFunctionMatrix(t -> [0 1; -10*cos(t) -24-10*sin(t)],2pi);

julia> ce = psceig(A; reltol = 1.e-10)
2-element Vector{Float64}:
 -6.398432404426897
 -1.4291292368715818e-13
````

However, full accuracy can be achieved with multiple shooting by determining the monodromy matrix as a product of, say 500 matrices, and computing the eigenvalues using the periodic Schur decomposition:

````JULIA
julia> ce = psceig(A,500; reltol = 1.e-10)
2-element Vector{Float64}:
   3.180554681463513e-16
 -23.99999999998858
````

Note that the evaluation of the 500 factors can be done in parallel, in which case a substantial speedup of computations can be achieved.

## References

[1] A. Varga. [A Periodic Systems Toolbox for Matlab](https://elib.dlr.de/12283/1/varga_ifac2005p1.pdf). Proc. of IFAC 2005 World Congress, Prague, Czech Republic, 2005.

[2] S. Bittanti and P. Colaneri. Periodic Systems - Filtering and Control, Springer Verlag, 2009.

[3] J. A. Richards. Analysis of Periodically Time-Varying Systems, Springer Verlag, 1983.
