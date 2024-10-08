# Periodic matrix utilities

* **[`pseig`](@ref)**   Characteristic multipliers of a periodic matrix.
* **[`psceig`](@ref)**   Characteristic exponents of a periodic matrix.
* **[`pseigsm`](@ref)**   Characteristic multipliers of a periodic symbolic matrix.
* **[`psceighr`](@ref)**   Characteristic exponents of a periodic matrix in Harmonic Array representation.
* **[`psceigfr`](@ref)**   Characteristic exponents of a periodic matrix in Fourier Function Matrix representation.
* **[`psceigsm`](@ref)**   Characteristic exponents of a periodic matrix in symbolic representation.
* **[`monodromy`](@ref)**  Monodromy matrix of a linear periodic time-varying system of ODE.
* **[`tvstm`](@ref)**  State transition matrix of a linear periodic time-varying system of ODE.
* **[`tvmeval`](@ref)**  Time response evaluation of a continuous-time periodic matrix. 
* **[`hreval`](@ref)**  Evaluation of a harmonic array for a numerical or symbolic time value. 
* **[`hrchop`](@ref)**  Removal of the negligible trailing terms of a harmonic representation. 
* **[`hrtrunc`](@ref)**  Truncation of a harmonic representation.  
* **[`pmaverage`](@ref)**  Evaluation of the time averaged matrix of a continuous-time periodic matrix. 


```@docs
pseig(at::PM, K::Int64; lifting, solver, reltol, abstol, dt) where {T, PM<:Union{HarmonicArray{:c, T}, PeriodicFunctionMatrix{:c, T}}}
pseig(A::PeriodicArray{:d, T}; fast) where T
pseig(A::PeriodicMatrix{:d, T}, k::Int64; fast) where T
psceig
pseigsm
psceighr
psceigfr
psceigsm
monodromy
tvstm
tvmeval
hreval
hrchop
hrtrunc
pmaverage
```
