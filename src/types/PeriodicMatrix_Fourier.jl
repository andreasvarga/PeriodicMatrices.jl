
struct FourierFunctionMatrix{Domain,T,X} <: AbstractPeriodicArray{Domain,T} 
   M::X
   period::Float64
   nperiod::Int
end 
# struct FourierFunctionMatrix{Domain,T} <: AbstractPeriodicArray{Domain,T} 
#    M::Fun
#    period::Float64
#    nperiod::Int
# end
# additional constructors
"""
     FourierFunctionMatrix(Afun, T) -> A::FourierFunctionMatrix

Continuous-time Fourier function matrix representation.

The Fourier function matrix object `A` of period `T` is built using
the Fourier series representation of a periodic matrix `Afun(t)` of subperiod `T′ = T/k`, 
where each entry of `Afun(t)` has the form

             p
      a_0 +  ∑ ( ac_i*cos(i*t*2*π/T′)+as_i*sin(i*2*π*t/T′) ) ,
            i=1 

where `k ≥ 1` is the number of subperiods (default: `k = 1`).   
The matrix `Afun` containing the Fourier representation, the period `T` and the 
number of subperiods `k` can be accessed via `A.M`, `A.period` and `A.nperiod`, respectively.
"""
function FourierFunctionMatrix{:c,T}(A::Fun, period::Real) where {T}
   period > 0 || error("period must be positive") 
   n, m = size(A,1), size(A,2)
   sint = domain(A)
   (sint.a == 0 && sint.b > 0) || error("the domain must be of the form 0..period")
   ti = rationalize(period/sint.b)
   denominator(ti) == 1 || error("only integer multiple periods supported")
   FourierFunctionMatrix{:c,eltype(domain(A)),Fun}(m == 1 ? reshape(A,n,m) : A, Float64(period), numerator(ti)) 
end
FourierFunctionMatrix(A::Fun, period::Real; nperiod::Int = 1)  = 
       FourierFunctionMatrix{:c,eltype(domain(A)),Fun}(A, period, nperiod) 
# FourierFunctionMatrix(A::Matrix{<:Fun}, period::Real; nperiod::Int = 1)  = 
#        FourierFunctionMatrix{:c,eltype(domain(A))}(A, period, nperiod) 

FourierFunctionMatrix(A::Fun)  = FourierFunctionMatrix{:c,eltype(domain(A))}(A::Fun, domain(A).b) 
function FourierFunctionMatrix{:c,T}(A::FourierFunctionMatrix, period::Real) where {T}
   period > 0 || error("period must be positive") 
   Aperiod = A.period
   r = rationalize(Aperiod/period)
   n, d = numerator(r), denominator(r)
   min(n,d) == 1 || error("new period is incommensurate with the old period")
   if period >= Aperiod
      FourierFunctionMatrix{:c,T,Fun}(A.M, Aperiod*d, A.nperiod*d)
   elseif period < Aperiod
      nperiod = div(A.nperiod,n)
      nperiod < 1 && error("new period is incommensurate with the old period")
      FourierFunctionMatrix{:c,T,Fun}(A.M, Aperiod/n, nperiod)
   end
end
set_period(A::FourierFunctionMatrix, period::Real) = FourierFunctionMatrix{:c,eltype(A)}(A,period)

FourierFunctionMatrix{:c,T}(A0::VecOrMat, period::Real) where {T <: Real}  = 
    FourierFunctionMatrix{:c,float(T)}(Fun(t->float(T).(A0),Fourier(0..period)), period) 
FourierFunctionMatrix(A0::VecOrMat{T}, period::Real) where {T <: Real}  = 
    FourierFunctionMatrix{:c,float(T)}(A0, period) 
function isconstant(A::FourierFunctionMatrix)
   for i = 1:size(A,1)
       for j = 1: size(A,2)
           ncoefficients(chop(A.M[i,j])) <= 1 || (return false)
       end
   end
   return true
end
# function isconstant(A::FourierFunctionMatrix)
#    k = length(size(A.M))
#    if  k > 1
#      for i = 1:size(A,1)
#           for j = 1: size(A,2)
#               ncoefficients(chop(A.M[i,j])) <= 1 || (return false)
#           end
#       end
#    elseif k == 1
#       for i = 1:size(A,1)
#           ncoefficients(chop(A.M[i])) <= 1 || (return false)
#       end
#    else
#       ncoefficients(chop(A.M)) <= 1 || (return false)      
#    end
#    return true
# end

Base.size(A::FourierFunctionMatrix) = (size(A.M,1),size(A.M,2))
Base.size(A::FourierFunctionMatrix, d::Integer) = d <= 2 ? size(A)[d] : 1
Base.eltype(A::FourierFunctionMatrix{:c,T}) where T = T

function Base.getindex(A::PM, inds...) where PM <: FourierFunctionMatrix
   size(inds, 1) != 2 &&
       error("Must specify 2 indices to index a periodic matrix")
   rows, cols = index2range(inds...) 
   FourierFunctionMatrix{:c,eltype(A),Fun}(A.M[rows,cols], A.period, A.nperiod)
end
function Base.lastindex(A::PM, dim::Int) where PM <: FourierFunctionMatrix
   size(A.M,dim)
end

# conversion to periodic function matrix
function Base.convert(::Type{PeriodicFunctionMatrix{:c,T}}, A::FourierFunctionMatrix) where T
   return PeriodicFunctionMatrix{:c,T}(x -> T.(A.M(T(x))), A.period, size(A), A.nperiod, isconstant(A))
end
Base.convert(::Type{PeriodicFunctionMatrix}, A::FourierFunctionMatrix) = convert(PeriodicFunctionMatrix{:c,eltype(A)}, A)

# function Base.convert(::Type{PeriodicFunctionMatrix}, A::FourierFunctionMatrix)
#    return PeriodicFunctionMatrix{:c,eltype(A)}(x -> A.M(x), A.period, size(A.M), A.nperiod, isconstant(A))
# end

# conversions to continuous-time Fourier function matrix
function Base.convert(::Type{FourierFunctionMatrix}, A::PeriodicFunctionMatrix) 
   T = float(eltype(A))
   return FourierFunctionMatrix{:c,T,Fun}(Fun(x -> T.(A.f(x)), Fourier(0..A.period/A.nperiod)), Float64(A.period), A.nperiod)
end
function Base.convert(::Type{FourierFunctionMatrix}, A::HarmonicArray) 
   tA = convert(PeriodicFunctionMatrix,A)
   return FourierFunctionMatrix{:c,eltype(tA),Fun}(Fun(x -> tA.f(x), Fourier(0..tA.period/tA.nperiod)), Float64(tA.period), tA.nperiod)
end
function Base.convert(::Type{FourierFunctionMatrix}, A::PeriodicTimeSeriesMatrix) 
   tA = convert(PeriodicFunctionMatrix,A)
   return FourierFunctionMatrix{:c,eltype(tA),Fun}(Fun(x -> tA.f(x), Fourier(0..tA.period/tA.nperiod)), Float64(tA.period), tA.nperiod)
end

# conversion to continuous-time HarmonicArray
Base.convert(::Type{HarmonicArray}, A::FourierFunctionMatrix) = ffm2hr(A)

# conversions to PeriodicTimeSeriesMatrix
function Base.convert(::Type{PeriodicTimeSeriesMatrix}, A::FourierFunctionMatrix; ns::Int = 128)
   ns > 0 || throw(ArgumentError("number of samples must be positive, got $ns"))
   ts = (0:ns-1)*A.period/ns/A.nperiod
   PeriodicTimeSeriesMatrix(tpmeval.(Ref(A),ts), A.period; nperiod = A.nperiod)
   #convert(PeriodicTimeSeriesMatrix,convert(PeriodicFunctionMatrix,A);ns)
end
