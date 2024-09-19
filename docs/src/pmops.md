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
horzcat
vertcat
blockdiag
iszero
isconstant
Base.isequal
Base.isapprox
issymmetric
```

## Operations on continuous-time periodic matrices


```@docs

LinearAlgebra.norm(A::PeriodicFunctionMatrix, p::Real; rtol)
trace(A::Union{HarmonicArray, PeriodicFunctionMatrix}; rtol)
pmderiv
```


## Operations on discrete-time periodic matrices


```@docs
Base.reverse
LinearAlgebra.norm(A::PeriodicArray, p::Real)
trace(A::PeriodicArray)
pmshift
pmsymadd!
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
