# Periodic matrix conversions

* **[`convert`](@ref)**   Conversion between discrete-time and between continuous-time periodic matrix representations.
* **[`ts2hr`](@ref)**   Conversion of  a periodic time series matrix to a harmonic array approximation.
* **[`pfm2hr`](@ref)**  Conversion of  a periodic function matrix to a harmonic array representation. 
* **[`ts2pfm`](@ref)**  Conversion of  an interpolated periodic time series matrix to a periodic function matrix.
* **[`hr2psm`](@ref)**  Conversion of  a harmonic array representation to a symbolic matrix.
* **[`psm2hr`](@ref)**  Conversion of  a periodic symbolic matrix into a harmonic array representation.
* **[`pm2pa`](@ref)**   Conversion of  a discrete-time periodic matrix object to a periodic array object.
* **[`ffm2hr`](@ref)**  Conversion of  a Fourier function matrix to a harmonic array representation. 
* **[`ffm2psm`](@ref)**  Conversion of a Fourier function matrix to a symbolic matrix. 
* **[`hr2bt`](@ref)**   Building a block Toeplitz matrix approximation of a harmonic (Fourier) array representation. 
* **[`hr2btupd`](@ref)**  Building an updated block Toeplitz matrix approximation of a harmonic (Fourier) array representation. 


```@docs
convert(::Type{<:PeriodicMatrix}, A::PeriodicArray{:d, T}) where T
convert(::Type{PeriodicFunctionMatrix}, ahr::HarmonicArray)
ts2hr
pfm2hr
ts2pfm
hr2psm
psm2hr
pm2pa
ffm2hr
ffm2psm
hr2bt
hr2btupd
```
