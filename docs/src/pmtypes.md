# Constructors for periodic matrices

## Discrete-time periodic matrices

* **[`PeriodicMatrix`](@ref)**   Discrete-time periodic matrix representation.
* **[`PeriodicArray`](@ref)**    Discrete-time periodic array representation.
* **[`SwitchingPeriodicMatrix`](@ref)** Discrete-time switching periodic matrix representation.
* **[`SwitchingPeriodicArray`](@ref)** Discrete-time switching periodic array representation.


```@docs
PeriodicMatrix
SwitchingPeriodicMatrix
size(A::PeriodicMatrix)
length(A::PeriodicMatrix)
eltype(A::PeriodicMatrix)
getindex(A::PeriodicMatrix, ind::Int64)
getindex(A::PeriodicMatrix, inds...)
lastindex(A::PeriodicMatrix)
lastindex(A::PeriodicMatrix, dim::Int64)
PeriodicArray
SwitchingPeriodicArray
size(A::PeriodicArray)
length(A::PeriodicArray)
eltype(A::PeriodicArray) 
getindex(A::PM, ind::Int64) where PM<:PeriodicArray
getindex(A::PeriodicArray, inds...)
lastindex(A::PeriodicArray)
lastindex(A::PeriodicArray, dim::Int64)
```

## Continuous-time periodic matrices

* **[`PeriodicFunctionMatrix`](@ref)**  Continuous-time periodic function matrix representation.
* **[`PeriodicSymbolicMatrix`](@ref)**   Continuous-time periodic symbolic matrix representation.
* **[`PeriodicTimeSeriesMatrix`](@ref)**   Continuous-time periodic time series matrix representation.
* **[`HarmonicArray`](@ref)**   Continuous-time harmonic array representation.
* **[`FourierFunctionMatrix`](@ref)**   Continuous-time Fourier functin matrix representation.
* **[`PeriodicSwitchingMatrix`](@ref)** Continuous-time switching periodic matrix representation.

```@docs
PeriodicFunctionMatrix
PeriodicSymbolicMatrix
PeriodicTimeSeriesMatrix
HarmonicArray
HarmonicArray(A0::MT, Acos::Union{Nothing, Vector{MT}}, Asin::Union{Nothing, Vector{MT}}, period::Real) where {T<:Real, MT<:VecOrMat{T}} 
FourierFunctionMatrix
PeriodicSwitchingMatrix
size(A::PeriodicFunctionMatrix)
eltype(A::PeriodicFunctionMatrix) 
getindex(A::PeriodicFunctionMatrix, inds...)
lastindex(A::PeriodicFunctionMatrix, dim::Int64)
```
