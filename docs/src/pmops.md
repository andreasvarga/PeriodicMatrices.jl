## Operations on periodic matrices

```@docs
Base.:-(A::PeriodicArray)
LinearAlgebra.inv
Base.transpose
Base.adjoint
opnorm
tr
Base.:+
Base.:-(A::PeriodicArray, B::PeriodicArray)
Base.:*
pmcopy
horzcat
vertcat
blockdiag
blockut
iszero
isconstant
Base.isequal
Base.isapprox
issymmetric
```

## Operations on continuous-time periodic matrices


```@docs

norm(A::Union{HarmonicArray, PeriodicFunctionMatrix}, p::Real; rtol)
trace(A::Union{HarmonicArray, PeriodicFunctionMatrix}; rtol)
pmderiv
pmrand(n::Int64, m::Int64, period::Real; nh)
```


## Operations on discrete-time periodic matrices


```@docs
Base.reverse
norm(A::PeriodicArray)
trace(A::PeriodicArray)
pmshift
pmsymadd!
pmata
pmaat
pmrand(::Type{PM}, m::Vector{Int64}, n::Vector{Int64}, period::Real) where PM<:PeriodicMatrix
pmrand(::Type{T}, m::Vector{Int64}, n::Vector{Int64}, period::Real) where T
```

## Operations with symmetric periodic matrices

```@docs
pmmulsym
pmtrmulsym
pmmultrsym
pmmuladdsym
pmmultraddsym
pmmuladdtrsym
```
