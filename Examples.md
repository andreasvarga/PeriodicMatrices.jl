# Examples

## Example 1: Floquet-analysis of differential equations with periodic parameters

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

This is the Example of Fig 3.3 in [3]. Assume again the period is $T = \pi$ and let $\psi(t)$ be the periodic negative slope sawtooth function $\psi(t) = -2t/T+1$. 
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

Stability can be illustrated with the lossy Hill equation with $\zeta = 0.2$:  

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
julia> ce 
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
          -10cos(t) & -24-10sin(t) 
          \end{array} \right]
```         
which has the characteristic exponents equal to $0$ and $-24$.

Using the standard settings to compute the characteristic exponents,  one of the resulting characteristic exponents has no one exact digit 
even by imposing a high relative tolerance for solving the underlying differential equations:

````JULIA
julia> using PeriodicMatrices

julia> A = PeriodicFunctionMatrix(t -> [0 1; -10*cos(t) -24-10*sin(t)],2pi);

julia> ce = psceig(A; reltol = 1.e-10)
2-element Vector{Float64}:
 -6.398432404426897
 -1.4291292368715818e-13
````

However, full accuracy can be achieved with the multiple shooting approach by determining the monodromy matrix as a product of, say 500, matrices and computing the 
characteristic multipliers (i.e., the eigenvalues of the monodromy matrix) using the _periodic Schur decomposition_:

````JULIA
julia> ce = psceig(A,500; reltol = 1.e-10)
2-element Vector{Float64}:
   3.180554681463513e-16
 -23.99999999998858
````

Note that the evaluation of the 500 factors can be done in parallel, in which case a substantial speedup of computations can be achieved.

Satisfactory accuracy can be also achieved using frequency lifting techniques based on the **FourierFunctionMatrix** representaion of $A(t)$:

````JULIA
julia> using PeriodicMatrices, ApproxFun

julia> A = FourierFunctionMatrix(Fun(t -> [0 1; -10*cos(t) -24-10*sin(t)], Fourier(0..2π)));

julia> ce = psceigfr(A,50)
2-element Vector{Float64}:
 -24.000000200907923
  -3.4668116010166546e-14
````
### Example 2: Discrete-time periodic matrices 

Let $A(t)$ be a continuous-time periodic matrix of period $T$, such that $A(t+T) = A(t)$ for all $t$. We can associate to $A(t)$ a discrete-time periodic matrix $A_d(k)$ of discrete period $p$ by defining
$A_d(k) = A((k-1)\Delta)$, where $\Delta := T/p$ is called the _sampling time_.  The sequence $A_d(1)$, ..., $A_d(k)$, ... satisfies $A_d(k+p) = A_d(k)$ for all $k$ and thus is $p$-periodic. 
When considering discrete-time periodic matrices, this connection to a continuous-time periodic matrix is implicitly assumed. 

A discrete-time periodic matrix $A_d(k)$ can be specified via a collection of component matrices  $A_1$, $A_2$, ..., $A_p$ (such that $A_k = A_d(k)$ ) and the real period $T$ (the discrete period $p$ implicitly results).
Normally, the component matrices have constant dimensions. However, for some specific problems, it is necessary
to allow for periodic time variability in the dimensions as well, in which case the component matrices $A_k$
exhibit a time-varying dimensionality.  
If the dimensions allow to form the product $\Phi(p,1) := A_p...A_2A_1$ such that $\Phi(p,1)$ is square, then $\Psi(1) := \Phi(p,1)$ is called the _monodromy matrix_ and, similarly to the continuous-time case,
its eigenvalues $\lambda_i$ are the _characteristic multipliers_ of $A_d(k)$.  The associated _characteristic exponents_ $\mu_i$ satisfy $\lambda_i = \mu_i^p$. 
The stability of a discrete-time periodic matrix can be assessed by computing its characteristic multipliers and checking that all characteristic multiplies have moduli less than one. 

The following examples are taken from [2], to illustrate some unexpected features. 
Consider first a discrete-time periodic matrix $A(k)$ of period $T = 2$ with two component matrices $A_1$ and $A_2$, with both having eigenvalues equal to zero. 

```math
   A_1 = \left[ \begin{array}{cc} 
          0 & 0\\
          2 & 0 
          \end{array} \right], \qquad 
   A_2 = \left[ \begin{array}{cc} 
          0 & 2\\
          0 & 0 
          \end{array} \right]

```         
In spite of this, the matrix is unstable as can be checked using a **PeriodicArray** representation of $A(k)$.


````JULIA
julia> using PeriodicMatrices 

julia> A = PeriodicArray(cat([0 0;2 0],[0 2;0 0], dims=3),2);

julia> pseig(A)
2-element Vector{Float64}:
 4.0
 0.0
````
Consider a discrete-time periodic matrix $A(k)$ of period $T = 2$ with two component matrices $A_1$ and $A_2$, with both having eigenvalues equal to 2. 

```math
   A_1 = \left[ \begin{array}{cc} 
          2 & 0\\
          -3.5 & 0 
          \end{array} \right], \qquad 
   A_2 = \left[ \begin{array}{cc} 
          2 & 1\\
          0 & 0 
          \end{array} \right]

```         
In spite of this, the matrix is stable as can be checked using a **PeriodicMatrix** representation of $A(k)$.


````JULIA
julia> using PeriodicMatrices

julia> A = PeriodicMatrix([[2 0;-3.5 0],[2 1;0 0]],2);

julia> pseig(A)
2-element Vector{Float64}:
 0.4999999999999999
 0.0
````
The above examples illustrate that there is no direct relation between the eigenvalues of component matrices and stability. 


Consider a discrete-time periodic matrix $A(k)$ of period $T = 2$ with time-varying dimensions with two component matrices $A_1$ and $A_2$ 

```math
   A_1 = \left[ \begin{array}{cc} 
          1 & 0
          \end{array} \right], \qquad 
   A_2 = \left[ \begin{array}{c} 
          1 \\
          1 
          \end{array} \right]

```         
The matrix has a _core_ characteristic multiplier equal to 1 and therefore can be considered marginally stable. 
The core characteristic multiplier can be isolated using the shifting of component matrices, which reverts the order of components used to form the monodromy matrix.


````JULIA
julia> using PeriodicMatrices

julia> A = PeriodicMatrix([[ 1 0 ], [ 1; 1]],2);

julia> pseig(A)
2-element Vector{Float64}:
 0.9999999999999997
 0.0

julia> pseig(pmshift(A,1))
1-element Vector{Float64}:
 0.9999999999999997

````


## References

[1] A. Varga. [A Periodic Systems Toolbox for Matlab](https://elib.dlr.de/12283/1/varga_ifac2005p1.pdf). Proc. of IFAC 2005 World Congress, Prague, Czech Republic, 2005.

[2] S. Bittanti and P. Colaneri. Periodic Systems - Filtering and Control, Springer Verlag, 2009.

[3] J. A. Richards. Analysis of Periodically Time-Varying Systems, Springer Verlag, 1983.
