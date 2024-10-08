# Utilities

* **[`pseig`](@ref)**   Eigenvalues of a cyclic matrix product.
* **[`ts2fm`](@ref)**   Compute the function matrix to interpolate a matrix time series.
* **[`psreduc_reg`](@ref)**  Fast reduction of a lifted regular pencil corresponding to a product of matrices. 

```@docs
pseig(A::Array{T, 3}; rev, fast) where T
pseig(A::Array{Matrix{T}, 1}, k::Int64; rev, fast) where T
ts2fm
psreduc_reg
```
