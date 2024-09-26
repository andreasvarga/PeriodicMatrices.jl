
"""
    PeriodicSymbolicMatrix(F, T; nperiod = k) -> A::PeriodicSymbolicMatrix

Continuous-time periodic symbolic matrix representation.
 
The continuous-time periodic symbolic matrix object `A` is built from `F`, a 
symbolic real matrix or vector of symbolic variable `t`, 
the associated time period `T` and the associated number of subperiods
specified via the keyword argument `nperiod = k`. 
It is assumed that  `F(t) = F(t+T/k)` for any real time value `t`.
The symbolic matrix `F`, the period `T` and the number of subperiods `k` 
can be accessed via `A.F`, `A.period` and `A.nperiod`, respectively.
"""
struct PeriodicSymbolicMatrix{Domain,T,X} <: AbstractPeriodicArray{Domain,T} 
   F::X
   period::Float64
   nperiod::Int
end 
# struct PeriodicSymbolicMatrix{Domain,T} <: AbstractPeriodicArray{Domain,T} 
#    F::Matrix{<:Num}
#    period::Float64
#    nperiod::Int
# end 

# additional constructors
function  PeriodicSymbolicMatrix{:c,T}(F::VecOrMat{T}, period::Real; nperiod::Int = 1) where {T <: Num} 
   period > 0 || error("period must be positive")       
   nperiod > 0 || error("number of subperiods must be positive") 
   # check that array F is depending only on t
   tt = rand()
   @variables t
   Ft = substitute.(F, (Dict(t => tt),))
   m, n = size(Ft,1), size(Ft,2)
   any(length.(Symbolics.get_variables.(Ft)) .> 0 ) && error("t must be the only variable in F")
   PeriodicSymbolicMatrix{:c,T,Matrix{T}}(n == 1 ? reshape(F,m,n) : F, Float64(period), nperiod) 
end
PeriodicSymbolicMatrix(F::VecOrMat{T}, period::Real; nperiod::Int = 1) where {T <: Union{Num,Real}} = 
    PeriodicSymbolicMatrix{:c,Num}(Num.(F), period; nperiod)
# PeriodicSymbolicMatrix(F::VecOrMat{T}, period::Real; nperiod::Int = 1) where {T <: Real} = 
#     PeriodicSymbolicMatrix{:c,Num}(Num.(F), Float64(period), nperiod)
# PeriodicSymbolicMatrix{:c,Num}(F::VecOrMat{T}, period::Real; nperiod::Int = 1) where {T <: Number} = 
#     PeriodicSymbolicMatrix{:c,Num}(Num.(F), Float64(period), nperiod)
function PeriodicSymbolicMatrix{:c,T}(A::PeriodicSymbolicMatrix, period::Real) where {T}
   period > 0 || error("period must be positive") 
   Aperiod = A.period
   r = rationalize(Aperiod/period)
   n, d = numerator(r), denominator(r)
   min(n,d) == 1 || error("new period is incommensurate with the old period")
   if period >= Aperiod
      PeriodicSymbolicMatrix{:c,T,Matrix{Num}}(A.F, Aperiod*d, A.nperiod*d)
   elseif period < Aperiod
      nperiod = div(A.nperiod,n)
      nperiod < 1 && error("new period is incommensurate with the old period")
      PeriodicSymbolicMatrix{:c,T,Matrix{Num}}(A.F, Aperiod/n, nperiod)
   end
end
set_period(A::PeriodicSymbolicMatrix, period::Real) = PeriodicSymbolicMatrix{:c,eltype(A)}(A,period)

# properties 
isconstant(A::PeriodicSymbolicMatrix) = all(length.(Symbolics.get_variables.(A.F)) .== 0)
# function isperiodic(A::VecOrMat{T}, period::Real) where {T <: Num} 
#    tt = rand()
#    @variables t
#    At = substitute.(A, (Dict(t => tt),))
#    return norm(At - substitute.(A, (Dict(t => tt+period),))) <= eps(10.)*max(1.,norm(At)) 
# end
#isperiodic(A::PeriodicSymbolicMatrix) = isconstant(A) ? true : isperiodic(A.F,A.period)
Base.size(A::PeriodicSymbolicMatrix) = size(A.F)
Base.size(A::PeriodicSymbolicMatrix, d::Integer) = d <= 2 ? size(A)[d] : 1
Base.eltype(A::PeriodicSymbolicMatrix{:c,T}) where T = T

function Base.getindex(A::PM, inds...) where PM <: PeriodicSymbolicMatrix
   size(inds, 1) != 2 &&
       error("Must specify 2 indices to index a periodic matrix")
   rows, cols = index2range(inds...) 
   PeriodicSymbolicMatrix{:c,eltype(A)}(A.F[rows,cols], A.period; nperiod = A.nperiod)
end
function Base.lastindex(A::PM, dim::Int) where PM <: PeriodicSymbolicMatrix
   lastindex(A.F,dim)
end


function Base.convert(::Type{PeriodicFunctionMatrix{:c,T}}, A::PeriodicSymbolicMatrix) where T
   @variables t
   f = eval(build_function(A.F, t, expression=Val{false}, nanmath=false)[1])
   PeriodicFunctionMatrix{:c,T}(x -> f(T(x)), A.period, size(A), A.nperiod, isconstant(A))
end
# function Base.convert(::Type{PeriodicFunctionMatrix}, A::PeriodicSymbolicMatrix) 
#    @variables t
#    f = eval(build_function(A.F, t, expression=Val{false}, nanmath=false)[1])
#    return PeriodicFunctionMatrix{:c,Float64}(x -> f(x), A.period, size(A), A.nperiod, isconstant(A))
# end
Base.convert(::Type{PeriodicFunctionMatrix}, A::PeriodicSymbolicMatrix) = convert(PeriodicFunctionMatrix{:c,Float64}, A)



# conversions to continuous-time PeriodicSymbolicMatrix
# function Base.convert(::Type{PeriodicSymbolicMatrix}, A::PeriodicFunctionMatrix) 
#    @variables t
#    # PeriodicSymbolicMatrix(Num.(A.f(t)), A.period; nperiod = A.nperiod)
#    M = try 
#       Num.(A.f(t))
#    catch 
#       return convert(PeriodicSymbolicMatrix,convert(HarmonicArray,A))
#    end
#    PeriodicSymbolicMatrix(M, A.period; nperiod = A.nperiod)
# end
function Base.convert(::Type{PeriodicSymbolicMatrix{:c,T}}, A::PeriodicFunctionMatrix) where {T}
   try 
      @variables t
      PeriodicSymbolicMatrix(Num.(A.f(t)), A.period; nperiod = A.nperiod)
   catch
      convert(PeriodicSymbolicMatrix,convert(HarmonicArray,A))
   end
end
Base.convert(::Type{PeriodicSymbolicMatrix}, A::PeriodicFunctionMatrix) = convert(PeriodicSymbolicMatrix{:c,eltype(A)}, A)
Base.convert(::Type{PeriodicSymbolicMatrix}, ahr::HarmonicArray) = 
   PeriodicSymbolicMatrix(hr2psm(ahr), ahr.period; nperiod = 1) 
   # PeriodicSymbolicMatrix(hr2psm(ahr), ahr.period; nperiod = ahr.nperiod)
Base.convert(::Type{PeriodicSymbolicMatrix}, A::PeriodicTimeSeriesMatrix) = 
   PeriodicSymbolicMatrix(hr2psm(convert(HarmonicArray,A)), A.period; nperiod =1)
   # PeriodicSymbolicMatrix(hr2psm(convert(HarmonicArray,A)), A.period; nperiod = A.nperiod)


# conversion to continuous-time HarmonicArray
Base.convert(::Type{HarmonicArray}, A::PeriodicSymbolicMatrix) = psm2hr(A)

# conversions to PeriodicTimeSeriesMatrix
# function Base.convert(::Type{PeriodicTimeSeriesMatrix}, A::PeriodicSymbolicMatrix; ns::Int = 128)
#    ns > 0 || throw(ArgumentError("number of samples must be positive, got $ns"))
#    tA = convert(PeriodicFunctionMatrix,A)
#    PeriodicTimeSeriesMatrix(tA.f.((0:ns-1)*tA.period/ns/tA.nperiod), tA.period; nperiod = tA.nperiod)
# end

function Base.convert(::Type{PeriodicTimeSeriesMatrix}, A::PeriodicSymbolicMatrix; ns::Int = 128)
   ns > 0 || throw(ArgumentError("number of samples must be positive, got $ns"))
   ts = (0:ns-1)*A.period/ns/A.nperiod
   PeriodicTimeSeriesMatrix(tpmeval.(Ref(A),ts), A.period; nperiod = A.nperiod)
end


function PeriodicFunctionMatrix(A::VecOrMat{Num}, period::Real) 
   @variables t
   f = eval(build_function(reshape(A,size(A,1),size(A,2)), t, expression=Val{false})[1])
   PeriodicFunctionMatrix{:c,Float64}(t -> f(t), period; isconst = false)
end

