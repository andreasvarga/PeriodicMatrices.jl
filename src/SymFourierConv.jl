
# conversions to continuous-time Fourier function matrix
function Base.convert(::Type{FourierFunctionMatrix}, A::PeriodicSymbolicMatrix) 
    # tA = convert(PeriodicFunctionMatrix,A)
    # return FourierFunctionMatrix{:c,eltype(tA),Fun}(Fun(x -> tA.f(x), Fourier(0..tA.period/tA.nperiod)), Float64(tA.period), tA.nperiod)
    return FourierFunctionMatrix{:c,Float64,Fun}(Fun(x -> tpmeval(A,x), Fourier(0..A.period/A.nperiod)), Float64(A.period), A.nperiod)
end  

function Base.convert(::Type{PeriodicSymbolicMatrix}, A::FourierFunctionMatrix)
    #convert(PeriodicSymbolicMatrix,convert(HarmonicArray,A))
    return PeriodicSymbolicMatrix{:c,Float64,Matrix{Num}}(ffm2psm(A), Float64(A.period), A.nperiod)
end
"""
     ffm2psm(Af::FourierFunctionMatrix, nrange atol = 0, rtol = 10*n*ϵ,) -> A::Matrix{Num}

Convert a range of harmonic components `nrange` of the Fourier function matrix `Af` to a symbolic matrix `A`. 
The default range is `nrange = 0:n`, where `n` is the order of the maximum harmonics. 
The tolerance used to assess nonzero coefficients is `tol = max(atol, rtol*maxnorm)`, where 
`maxnorm` is the maximum absolute value of the coefficients `ac_i` and `as_i` in `Af(t)`. The default values of tolerances
are `atol = 0` and `rtol = 10*n*ϵ`, where `ϵ` is the working machine precision.

"""
function ffm2psm(Af::FourierFunctionMatrix, nrange::Union{UnitRange,Missing} = missing; atol::Real = 0, rtol::Real = 0)  
   T = eltype(Af) 
   Symbolics.@variables t   
   Period = Af.period
   lens = length.(coefficients.(Matrix(Af.M)))
   n = max(div(maximum(lens)-1,2),0)+1
   n1, n2 = size(Af)
   if ismissing(nrange)
      i1 = 0; i2 = n
   else
      i1 = max(first(nrange),0)
      i2 = min(last(nrange), n)
      i1 = min(i1,i2)
   end
   a = zeros(Num,n1,n2)
   
   ts = t*2*pi*Af.nperiod/Period
   tol = iszero(atol) ? (iszero(rtol) ? 10*n*eps(maximum(norm.(coefficients.(Matrix(Af.M)),Inf))) : rtol*maximum(norm.(coefficients.(Matrix(Af.M)),Inf)) ) : atol


   for i = 1:n1
      for j = 1:n2
          tt = coefficients(Af.M[i,j])
          tt[abs.(tt) .<= tol] .= zero(T)
          lentt = 0
          for ii = lens[i,j]:-1:1
              iszero(tt[ii]) || (lentt = ii; break)
          end
          i1 == 0 && (a[i,j] = tt[1]) 
          if min(lentt,i2) > 0
             i1 == 0 ? k = 0 : k = i1-1
             for k1 = 2*max(i1,1):2:lentt
                 k += 1
                 a[i,j] += tt[k1]*(sin(k*ts))
             end
             i1 == 0 ? k = 0 : k = i1-1
             for k1 = 2*max(i1,1):2:lentt-1
                 k += 1
                 a[i,j] += tt[k1+1]*(cos(k*ts)) 
             end
          end
      end
   end
   return a
end 
