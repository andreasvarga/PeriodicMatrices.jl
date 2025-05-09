bldiag(mats...) = cat(mats..., dims=(1,2))
function promote_period(PM1,args...; ndigits = 4)
    nlim = 2^ndigits
    if typeof(PM1) <: AbstractVecOrMat 
       period = nothing 
    else
       period = PM1.period
    end
    # determine period
    for a in args
        typeof(a) <: AbstractVecOrMat && continue
        if isnothing(period) 
            period = a.period
           continue
        else
           peri = a.period
           r = rationalize(period/peri)
           num = numerator(r)
           den = denominator(r)
           num <= nlim || den <= nlim || error("incommensurate periods")
           period = period*den
        end
    end
    return period
end
function promote_period2(PM1,args...; ndigits = 4)
    #determine period
    period = promote_period(PM1,args...; ndigits)   
    isnothing(period) && (return nothing, nothing)
    # determine nperiod
    if typeof(PM1) <: AbstractVecOrMat 
       nperiod = nothing 
    else
       nperiod = PM1.nperiod*rationalize(period/PM1.period).num
    end
    for a in args
        typeof(a) <: AbstractVecOrMat && continue
        if isnothing(nperiod) 
           nperiod = a.nperiod*rationalize(period/a.period).num
        else
           nperiod = gcd(nperiod,a.nperiod*rationalize(period/a.period).num)
        end
    end
    return period, nperiod
end
 
# Operations with periodic arrays
"""
    pmshift(A[,k = 1])

Circular shifting of the components of a discrete-time periodic matrix `A` with `k` positions. 
"""
function pmshift(A::PeriodicArray, k::Int = 1)
    return PeriodicArray(A.M[:,:,mod.(k:A.dperiod+k-1,A.dperiod).+1], A.period; nperiod = A.nperiod)
end
# Inverse of a periodic matrices
"""
    inv(A)

Periodic matrix inversion.  Computes the periodic matrix `Ainv` such that `A * Ainv = Ainv*A = I`, where `I` is the identity matrix. 
"""
function LinearAlgebra.inv(A::PeriodicArray) 
    x = similar(A.M)
    [x[:,:,i] = inv(A.M[:,:,i]) for i in 1:size(A.M,3)]
    return PeriodicArray(x, A.period; nperiod = A.nperiod)
end
# Transpose of a periodic matrices
"""
    transpose(A)

Transpose of a periodic matrix.  
"""
function LinearAlgebra.transpose(A::PeriodicArray)
    return PeriodicArray(permutedims(A.M,(2,1,3)), A.period; nperiod = A.nperiod)
end
"""
    adjoint(A) or A'

Adjoint of a periodic matrix (equivalent to transpose). 
"""
function LinearAlgebra.adjoint(A::PeriodicArray)
    return PeriodicArray(permutedims(A.M,(2,1,3)), A.period; nperiod = A.nperiod)
end
"""
    reverse(A) 

Reverse the order of elements of a discrete-time periodic matrix. 
"""
function Base.reverse(A::PeriodicArray)
    return PeriodicArray(reverse(A.M,dims=3), A.period; nperiod = A.nperiod)
end
"""
    opnorm(A,p::Real=2) -> Anorm

Compute the point-wise time-dependent operator norm (or matrix norm) induced by the vector `p`-norm, where valid values of `p` are `1`, `2`, or `Inf`. 
For a periodic matrix `A(t)`, the resulting `Anorm(t)` is a periodic vector of same type as `A`, such that `Anorm(t)` is the `p`-norm of `A(t)` for each `t`. 
"""
function LinearAlgebra.opnorm(A::PeriodicArray, p::Union{Real,Missing} = missing)
    k = size(A.M,3)
    x = Array{eltype(A),3}(undef, 1, 1, k)
    if ismissing(p)
        [x[1,1,i] = norm(view(A.M,:,:,i)) for i in 1:k]  # Frobenius noorm
    else
        [x[1,1,i] = opnorm(view(A.M,:,:,i),p) for i in 1:k] # p-norm
    end
    return PeriodicArray(x,A.period; nperiod = A.nperiod)
end
"""
    norm(A,p::Real=2) -> Anorm

Compute the Frobenius-norm of a discrete-time periodic matrix. 
For a discrete-time periodic matrix `A(t)`, the resulting `Anorm` is the `p`-norm of the vector of Frobenius norms of each component of `A`. 
"""
function LinearAlgebra.norm(A::PeriodicArray, p::Real = 2)
    nrm = norm([norm(view(A.M,:,:,i)) for i in 1:size(A.M,3)],p)
    if p == 2 
       return nrm*sqrt(A.nperiod)
    elseif isinf(p) 
       return nrm
    elseif p == 1 
       return nrm*A.nperiod
    else 
        throw(ArgumentError("only p-norms for p = 1, 2, or Inf are supported"))
    end
end
"""
    tr(A) -> Atr

Compute the point-wise trace of periodic matrix. 
For a periodic matrix `A(t)`, the resulting `Atr(t)` is a periodic vector of same type as `A`, such that `Atr(t)` is the trace of `A(t)` for each `t`. 
"""
function LinearAlgebra.tr(A::PeriodicArray)
    p = size(A.M,3)
    x = Array{eltype(A),3}(undef, 1, 1, p)
    [x[1,1,i] = tr(view(A.M,:,:,i)) for i in 1:p]
    return PeriodicArray(x,A.period; nperiod = A.nperiod)
end
"""
    trace(A) -> Atrace

Compute the trace of periodic matrix. 
For a periodic matrix `A(t)`, the resulting `Atrace` is the sum of point-wise traces over a complete period (see [`tr`](@ref). 
"""
function trace(A::PeriodicArray)
    t = zero(eltype(A.M))
    for i in 1:size(A.M,3)
        t += tr(view(A.M,:,:,i))
    end
    return t*A.nperiod
end
# Addition of periodic matrices
"""
    +(A, B)

Periodic matrix addition. One of arguments may be a constant matrix or a uniform scaling. 
`A` and `B` may have different, but commensurate periods.

When adding continuous-time periodic mtrices, if one of arguments is of type `PeriodicFunctionMatrix`, 
then the result is also of type `PeriodicFunctionMatrix`. 
"""
function +(A::PeriodicArray, B::PeriodicArray)
    isconstant(A) || isconstant(B) || A.Ts ≈ B.Ts || error("A and B must have the same sampling time")
    period = promote_period(A, B)
    m, n, pa = size(A.M)
    mb, nb, pb = size(B.M)
    (m, n) == (mb, nb) || throw(DimensionMismatch("A and B must have the same dimensions"))
    nta = numerator(rationalize(period/A.period))
    ntb = numerator(rationalize(period/B.period))
    pa == pb && nta == 1 && ntb == 1 && (return PeriodicArray(A.M+B.M, A.period; nperiod = A.nperiod))
    K = max(nta*A.nperiod*pa,ntb*B.nperiod*pb)
    #K = nta*A.nperiod*pa
    p = lcm(pa,pb)
    T = promote_type(eltype(A),eltype(B))
    X = Array{T,3}(undef, m, n, p)
    for i = 1:p
        ia = mod(i-1,pa)+1
        ib = mod(i-1,pb)+1
        X[:,:,i] = A.M[:,:,ia]+B.M[:,:,ib]
    end
    return PeriodicArray(X, period; nperiod = div(K,p))
end
+(A::PeriodicArray, C::AbstractMatrix) = +(A, PeriodicArray(C, A.Ts; nperiod = 1))
+(A::AbstractMatrix, C::PeriodicArray) = +(PeriodicArray(A, C.Ts; nperiod = 1), C)
# Unary minus operator
"""
    -(A)

Unary minus operator. 
"""
function -(A::PeriodicArray)
    PeriodicArray(-A.M, A.period; nperiod = A.nperiod)
end
# Substracion of periodic matrices
"""
    -(A, B)

Periodic matrix substraction. One of arguments may be a constant matrix or a uniform scaling. 
`A` and `B` may have different, but commensurate periods.

When substracting continuous-time periodic mtrices, if one of arguments is of type `PeriodicFunctionMatrix`, 
then the result is also of type `PeriodicFunctionMatrix`. 
"""
function -(A::PeriodicArray, B::PeriodicArray)
    +(A,-B)
end
-(A::PeriodicArray, C::AbstractMatrix) = +(A,-C)
-(A::AbstractMatrix, C::PeriodicArray) = +(A, -C)
function (+)(A::PeriodicArray, J::UniformScaling{<:Real})
    m, n = size(A,1), size(A,2)
    n == m || throw(DimensionMismatch("periodic matrix is not square: dimensions are $((m,n))"))
    x = similar(A.M)
    [x[:,:,i] = A.M[:,:,i] + J for i in 1:size(A.M,3)]
    return PeriodicArray(x, A.period; nperiod = A.nperiod)
end
(+)(J::UniformScaling{<:Real}, A::PeriodicArray) = +(A,J)
(-)(A::PeriodicArray, J::UniformScaling{<:Real}) = +(A,-J)
(-)(J::UniformScaling{<:Real}, A::PeriodicArray) = +(-A,J)
"""
    pmsymadd!(A, α = 1)

Compute for a discrete-time periodic matrix `A` the symmetric matrix `α*(A+A')` in place. 
"""
function pmsymadd!(A::Union{PeriodicArray,SwitchingPeriodicArray}, α = 1)
    # compute the symmetric matrix α*(A+transpose(A))
    m, n = size(A) 
    m == n || throw(ArgumentError("matrix A must be square"))
    if α == 1
       for i = 1:length(A) 
           inplace_transpose_add!(view(A.M,:,:,i))
       end
    else
        for i = 1:length(A) 
            inplace_transpose_add!(view(A.M,:,:,i),α)
        end
    end
    return A
end

# Multiplication of periodic matrices
"""
    *(A, B)

Periodic matrix multiplication. One of arguments may be a constant matrix, a scalar or a uniform scaling. 
`A` and `B` may have different, but commensurate periods.

Whenmultiplying continuous-time periodic mtrices, if one of arguments is of type `PeriodicFunctionMatrix`, 
then the result is also of type `PeriodicFunctionMatrix`. 
"""
function *(A::PeriodicArray, B::PeriodicArray)
    isconstant(A) || isconstant(B) || A.Ts ≈ B.Ts || error("A and B must have the same sampling time")
    period = promote_period(A, B)
    m, na, pa = size(A.M)
    mb, n, pb = size(B.M)
    na == mb || throw(DimensionMismatch("A and B have incompatible dimensions"))
    nta = numerator(rationalize(period/A.period))
    ntb = numerator(rationalize(period/B.period))
    K = max(nta*A.nperiod*pa,ntb*B.nperiod*pb)
    #K == B.nperiod*pb || error("A and B must have the same sampling time")
    #pa == pb && (return PeriodicArray(A.M*B.M, A.period; nperiod = A.nperiod))
    p = lcm(pa,pb)
    T = promote_type(eltype(A),eltype(B))
    X = Array{T,3}(undef, m, n, p)
    for i = 1:p
        ia = mod(i-1,pa)+1
        ib = mod(i-1,pb)+1
        mul!(view(X,:,:,i), view(A.M,:,:,ia), view(B.M,:,:,ib))
    end
    return PeriodicArray(X, period; nperiod = div(K,p))
end
*(A::PeriodicArray, C::AbstractMatrix) = *(A, PeriodicArray(C, A.Ts; nperiod = 1))
*(A::AbstractMatrix, C::PeriodicArray) = *(PeriodicArray(A, C.Ts; nperiod = 1), C)
*(A::PeriodicArray, C::Real) = PeriodicArray(C*A.M, A.period; nperiod = A.nperiod)
*(A::Real, C::PeriodicArray) = PeriodicArray(A*C.M, C.period; nperiod = C.nperiod)
/(A::PeriodicArray, C::Real) = *(A, 1/C)
*(J::UniformScaling{<:Real}, A::PeriodicArray) = J.λ*A
*(A::PeriodicArray, J::UniformScaling{<:Real}) = A*J.λ

for (PMF, MF) in ((:pmmuladdsym, :muladdsym!), (:pmmultraddsym, :multraddsym!), (:pmmuladdtrsym,:muladdtrsym!) )
    @eval begin
        function $PMF(A::PeriodicArray,B::PeriodicArray,C::PeriodicArray, (α,β) = (true, true))
            if isconstant(A)
               isconstant(B) || isconstant(C) || B.Ts ≈ C.Ts || error("B and C must have the same sampling time")
            elseif isconstant(B)
               isconstant(C) || A.Ts ≈ C.Ts || error("A and C must have the same sampling time")
            elseif isconstant(C)
                A.Ts ≈ B.Ts || error("A and B must have the same sampling time") 
            else
                A.Ts ≈ B.Ts ≈ C.Ts || error("A, B and C must have the same sampling time")
            end
            period = promote_period(A, B, C)
            pa = length(A) 
            pb = length(B)
            pc = length(C)
            p = lcm(pa,pb,pc)
            T = promote_type(eltype(A),eltype(B),eltype(C))
            n = size(A,1)
            X = Array{T,3}(undef, n, n, p)
            nta = numerator(rationalize(period/A.period))
            ntb = numerator(rationalize(period/B.period))
            ntc = numerator(rationalize(period/C.period))
            Ka = nta*A.nperiod*pa
            Kb = ntb*B.nperiod*pb
            Kc = ntc*C.nperiod*pc
            K = max(Ka,Kb,Kc)
        
            for i = 1:p
                ia = mod(i-1,pa)+1
                ib = mod(i-1,pb)+1
                ic = mod(i-1,pc)+1
                copyto!(view(X,:,:,i), view(A.M,:,:,ia))
                $MF(view(X,:,:,i), view(B.M,:,:,ib), view(C.M,:,:,ic),(α,β))
            end
            return PeriodicArray(X, period; nperiod = div(K,p))
        end
        $PMF(A::PeriodicArray, B::PeriodicArray, C::PeriodicArray, α, β) = $PMF(A, B, C, (α,β))
        function $PMF(A::AbstractMatrix,B::AbstractMatrix,C::AbstractMatrix, (α,β) = (true, true)) 
            T = promote_type(eltype(A),eltype(B),eltype(C))
            $MF(LinearAlgebra.copy_oftype(A,T),B,C,(α,β))
        end
        $PMF(A::AbstractMatrix,B::AbstractMatrix,C::AbstractMatrix, α, β) = $PMF(A, B, C, (α,β))
     end
end
"""
    pmmuladdsym(A, B, C, α, β)

Compute the symmetric periodic matrix `α*A + β*B*C`, where `α` and `β` are real scalars, 
`A` is a symmetrix periodic matrix and the product `B*C` is known to be symmetric. 
All matrix arguments may be constant matrices as well.

_Note:_ This function is available only for periodic matrices of types `PeriodicArray`, `PeriodicMatrix`, `PeriodicFunctionMatrix` and `HarmonicArray`.
"""
pmmuladdsym(A::PeriodicArray, B::PeriodicArray, C::PeriodicArray, α::Real, β::Real)
"""
    pmmultraddsym(A, B, C, α, β)

Compute the symmetric periodic matrix `α*A + β*B'*C`, where `α` and `β` are real scalars, 
`A` is a symmetrix periodic matrix and the product `B'*C` is known to be symmetric. 
All matrix arguments may be constant matrices as well.

_Note:_ This function is available only for periodic matrices of types `PeriodicArray`, `PeriodicMatrix`, `PeriodicFunctionMatrix` and `HarmonicArray`.
"""
pmmultraddsym(A::PeriodicArray, B::PeriodicArray, C::PeriodicArray, α::Real, β::Real)
"""
    pmmuladdtrsym(A, B, C, α, β)

Compute the symmetric periodic matrix `α*A + β*B*C'`, where `α` and `β` are real scalars, 
`A` is a symmetrix periodic matrix and the product `B*C'` is known to be symmetric. 
All matrix arguments may be constant matrices as well.

_Note:_ This function is available only for periodic matrices of types `PeriodicArray`, `PeriodicMatrix`, `PeriodicFunctionMatrix` and `HarmonicArray`.
"""
pmmuladdtrsym(A::PeriodicArray, B::PeriodicArray, C::PeriodicArray, α::Real, β::Real)


for (PMF, MF) in ((:pmmulsym, :muladdsym!), (:pmtrmulsym, :multraddsym!), (:pmmultrsym,:muladdtrsym!) )
    @eval begin
        function $PMF(B::PeriodicArray,C::PeriodicArray, β = true)
            # compute the symmetric results βB*C, βB'*C, and βB*C'. 
            isconstant(B) || isconstant(C) || B.Ts ≈ C.Ts || error("B and C must have the same sampling time")
            period = promote_period(B, C)
            pb = length(B)
            pc = length(C)
            p = lcm(pb,pc)
            T = promote_type(eltype(B),eltype(C))
            n = $PMF == pmtrmulsym ? size(B,2) : size(B,1)
            X = Array{T,3}(undef, n, n, p)
            ntb = numerator(rationalize(period/B.period))
            ntc = numerator(rationalize(period/C.period))
            K = max(ntb*B.nperiod*pb,ntc*C.nperiod*pc)   
            for i = 1:p
                ib = mod(i-1,pb)+1
                ic = mod(i-1,pc)+1
                $MF(view(X,:,:,i), view(B.M,:,:,ib), view(C.M,:,:,ic),(0,β))
            end
            return PeriodicArray(X, period; nperiod = div(K,p))
        end
    end
end
"""
    pmmulsym(B, C, β)

Compute the symmetric periodic matrix `β*B*C`, where `β` is a real scalar 
and the product `B*C` is known to be symmetric. 
All matrix arguments may be constant matrices as well.

_Note:_ This function is available only for periodic matrices of types `PeriodicArray`, `PeriodicMatrix`, `PeriodicFunctionMatrix` and `HarmonicArray`.
"""
pmmulsym(B::PeriodicArray, C::PeriodicArray, β::Real)
"""
    pmtrmulsym(B, C, β)

Compute the symmetric periodic matrix `β*B'*C`, where `β` is a real scalar
and the product `B'*C` is known to be symmetric. 
All matrix arguments may be constant matrices as well.

_Note:_ This function is available only for periodic matrices of types `PeriodicArray`, `PeriodicMatrix`, `PeriodicFunctionMatrix` and `HarmonicArray`.
"""
pmtrmulsym(B::PeriodicArray, C::PeriodicArray, β::Real)
"""
    pmmultrsym(B, C, β)

Compute the symmetric periodic matrix `β*B*C'`, where `β` is a real scalar 
and the product `B*C'` is known to be symmetric. 
All matrix arguments may be constant matrices as well.

_Note:_ This function is available only for periodic matrices of types `PeriodicArray`, `PeriodicMatrix`, `PeriodicFunctionMatrix` and `HarmonicArray`.
"""
pmmultrsym(B::PeriodicArray, C::PeriodicArray, β::Real)


for PMF in (:pmmuladdsym, :pmmultraddsym, :pmmuladdtrsym)
    for PM in (:PeriodicArray, :PeriodicMatrix)
        @eval begin
            $PMF(A::$PM,B::AbstractMatrix,C::$PM, (α,β) = (true, true)) = $PMF(A, $PM(B, A.period), C, (α,β))
            $PMF(A::$PM,B::AbstractMatrix,C::$PM, α, β) = $PMF(A, B, C, (α,β))
            $PMF(A::$PM,B::$PM,C::AbstractMatrix, (α,β) = (true, true)) = $PMF(A, B, $PM(C, A.period), (α,β))
            $PMF(A::$PM,B::$PM,C::AbstractMatrix, α, β) = $PMF(A, B, C, (α,β))
            $PMF(A::$PM,B::AbstractMatrix,C::AbstractMatrix, (α,β) = (true, true)) = $PMF(A, $PM(B, A.period), $PM(C, A.period), (α,β))
            $PMF(A::$PM,B::AbstractMatrix,C::AbstractMatrix, α, β) = $PMF(A, B, C, (α,β))
            $PMF(A::AbstractMatrix,B::$PM,C::$PM, (α,β) = (true, true)) = $PMF($PM(A, B.period), B, C, (α,β))
            $PMF(A::AbstractMatrix,B::$PM,C::$PM, α, β) = $PMF(A, B, C, (α,β))
            $PMF(A::AbstractMatrix,B::AbstractMatrix,C::$PM, (α,β) = (true, true)) = $PMF($PM(A, C.period), $PM(B, C.period), C, (α,β))
            $PMF(A::AbstractMatrix,B::AbstractMatrix,C::$PM, α, β) = $PMF(A, B, C, (α,β))
            $PMF(A::AbstractMatrix,B::$PM,C::AbstractMatrix, (α,β) = (true, true)) = $PMF($PM(A, B.period), B, $PM(C, B.period), (α,β))
            $PMF(A::AbstractMatrix,B::$PM,C::AbstractMatrix, α, β) = $PMF(A, B, C, (α,β))
        end
    end
end
for PMF in (:pmmulsym, :pmtrmulsym, :pmmultrsym)
    for PM in (:PeriodicArray, :PeriodicMatrix)
        @eval begin
            $PMF(B::$PM,C::AbstractMatrix, β = true) = $PMF(B, $PM(C, B.period), β)
            $PMF(B::AbstractMatrix,C::$PM, β = true) = $PMF($PM(B, C.period), C, β)
        end
    end
end
   
"""
    issymmetric(A)

Symmetry check of a periodic matrix.     
"""
function LinearAlgebra.issymmetric(A::PeriodicArray)
    all([issymmetric(A.M[:,:,i]) for i in 1:A.dperiod])
end
"""
    iszero(A)

Exact equality check with a zero periodic matrix.     
"""
function Base.iszero(A::PeriodicArray)
    iszero(A.M)
end
"""
    isequal(A, B)

Exact equality check of two periodic matries. `isequal(A,B)` is equivalent to the syntax `A == B`.
`A` and `B` must have equal subperiods.    
"""
function isequal(A::PeriodicArray, B::PeriodicArray)
    A.period == B.period && A.nperiod == B.nperiod && isequal(A.M, B.M)
end
function (==)(A::PeriodicArray, B::PeriodicArray)
    A.period == B.period && A.nperiod == B.nperiod && isequal(A.M, B.M)
end

"""
    isapprox(A, B; rtol::Real = sqrt(eps(Float64)), atol::Real = 0)

Inexact equality comparison of two periodic matries. `isaprox(A,B)` is equivalent to the syntax `A ≈ B`.
`A` and `B` must have equal subperiods. One of arguments may be an uniform scaling.   

`atol` and `rtol` are the absolute tolerance and relative tolerance, respectively, to be used for comparison.
"""
function Base.isapprox(A::PeriodicArray, B::PeriodicArray; rtol::Real = sqrt(eps(Float64)), atol::Real = 0)
    A.period == B.period && A.nperiod == B.nperiod && isapprox(A.M, B.M; rtol, atol)
end
function Base.isapprox(A::PeriodicArray, J::UniformScaling{<:Real}; kwargs...)
    all([isapprox(A.M[:,:,i], J; kwargs...) for i in 1:size(A.M,3)])
end
Base.isapprox(J::UniformScaling{<:Real}, A::PeriodicArray; kwargs...) = isapprox(A, J; kwargs...)


"""
    horzcat(A, B)
    hcat(A,B)

Horizontal concatenation of two periodic matrices. Equivalent to the syntax `[A B]`.
`A` and `B` may have different, but commensurate periods. One of arguments may be a constant matrix.
"""
function horzcat(A::PeriodicArray, B::PeriodicArray)
    A.Ts ≈ B.Ts || error("A and B must have the same sampling time")
    period = promote_period(A, B)
    m, na, pa = size(A.M)
    mb, nb, pb = size(B.M)
    m == mb || throw(DimensionMismatch("A and B have incompatible row dimensions"))
    p = lcm(pa,pb)
    nta = numerator(rationalize(period/A.period))
    ntb = numerator(rationalize(period/B.period))
    K = max(nta*A.nperiod*pa,ntb*B.nperiod*pb)
    T = promote_type(eltype(A),eltype(B))
    X = Array{T,3}(undef, m, na+nb, p)
    for i = 1:p
        ia = mod(i-1,pa)+1
        ib = mod(i-1,pb)+1
        X[:,:,i] = [view(A.M,:,:,ia) view(B.M,:,:,ib)]
    end
    return PeriodicArray(X, period; nperiod = div(K,p))
end
Base.hcat(A::PeriodicArray, B::PeriodicArray) = horzcat(A,B)
horzcat(A::PeriodicArray, B::AbstractMatrix) = horzcat(A, PeriodicArray(B, A.Ts; nperiod = 1))
horzcat(A::AbstractMatrix, B::PeriodicArray) = horzcat(PeriodicArray(A, B.Ts; nperiod = 1), B)
hcat(A::PeriodicArray, B::AbstractMatrix) = horzcat(A, PeriodicArray(B, A.Ts; nperiod = 1))
hcat(A::AbstractMatrix, B::PeriodicArray) = horzcat(PeriodicArray(A, B.Ts; nperiod = 1), B)

"""
    vertcat(A, B)
    vcat(A,B)

Vertical concatenation of two periodic matrices. Equivalent to the syntax `[A; B]`.
`A` and `B` may have different, but commensurate periods.
"""
function vertcat(A::PeriodicArray, B::PeriodicArray)
    A.Ts ≈ B.Ts || error("A and B must have the same sampling time")
    period = promote_period(A, B)
    ma, n, pa = size(A.M)
    mb, nb, pb = size(B.M)
    n == nb || throw(DimensionMismatch("A and B have incompatible column dimensions"))
    p = lcm(pa,pb)
    nta = numerator(rationalize(period/A.period))
    ntb = numerator(rationalize(period/B.period))
    K = max(nta*A.nperiod*pa,ntb*B.nperiod*pb)
    T = promote_type(eltype(A),eltype(B))
    X = Array{T,3}(undef, ma+mb, n, p)
    for i = 1:p
        ia = mod(i-1,pa)+1
        ib = mod(i-1,pb)+1
        X[:,:,i] = [view(A.M,:,:,ia); view(B.M,:,:,ib)]
    end
    return PeriodicArray(X, period; nperiod = div(K,p))
end
Base.vcat(A::PeriodicArray, B::PeriodicArray) = vertcat(A,B)
vertcat(A::PeriodicArray, B::AbstractMatrix) = vertcat(A, PeriodicArray(B, A.Ts; nperiod = 1))
vertcat(A::AbstractMatrix, B::PeriodicArray) = vertcat(PeriodicArray(A, B.Ts; nperiod = 1), B)
vcat(A::PeriodicArray, B::AbstractMatrix) = vertcat(A, PeriodicArray(B, A.Ts; nperiod = 1))
vcat(A::AbstractMatrix, B::PeriodicArray) = vertcat(PeriodicArray(A, B.Ts; nperiod = 1), B)

"""
    blockdiag(A, B)

Block diagonal appending of two periodic matrices. 
`A` and `B` may have different, but commensurate periods.
"""
function blockdiag(A::PeriodicArray, B::PeriodicArray)
    A.Ts ≈ B.Ts || error("A and B must have the same sampling time")
    period = promote_period(A, B)
    ma, na, pa = size(A.M)
    mb, nb, pb = size(B.M)
    p = lcm(pa,pb)
    nta = numerator(rationalize(period/A.period))
    K = nta*A.nperiod*pa
    T = promote_type(eltype(A),eltype(B))
    X = Array{T,3}(undef, ma+mb, na+nb, p)
    for i = 1:p
        ia = mod(i-1,pa)+1
        ib = mod(i-1,pb)+1
        X[:,:,i] = bldiag(view(A.M,:,:,ia), view(B.M,:,:,ib))
    end
    return PeriodicArray(X, period; nperiod = div(K,p))
end



# Operations with periodic matrices
function pmshift(A::PeriodicMatrix, k::Int = 1)
    return PeriodicMatrix(A.M[mod.(k:A.dperiod+k-1,A.dperiod).+1], A.period; nperiod = A.nperiod)
end
LinearAlgebra.inv(A::PeriodicMatrix) = PeriodicMatrix(inv.(A.M), A.period; nperiod = A.nperiod)
function LinearAlgebra.transpose(A::PeriodicMatrix)
    return PeriodicMatrix(copy.(transpose.(A.M)), A.period; nperiod = A.nperiod)
end
function LinearAlgebra.adjoint(A::PeriodicMatrix)
    return PeriodicMatrix(copy.(transpose.(A.M)), A.period; nperiod = A.nperiod)
end
function Base.reverse(A::PeriodicMatrix)
    return PeriodicMatrix(reverse(A.M), A.period; nperiod = A.nperiod)
end
function LinearAlgebra.tr(A::PeriodicMatrix)
    return PeriodicMatrix([[tr(A.M[i])] for i in 1:length(A)], A.period; nperiod = A.nperiod)
end
function trace(A::PeriodicMatrix)
    t = zero(eltype(A))
    for i in 1:length(A)
        t += tr(A.M[i])
    end
    return t*A.nperiod
end
function LinearAlgebra.opnorm(A::PeriodicMatrix, p::Union{Real,Missing} = missing)
    if ismissing(p)
        return PeriodicMatrix([norm(A.M[i]) for i in 1:length(A)], A.period; nperiod = A.nperiod)
    else
        return PeriodicMatrix([opnorm(A.M[i],p) for i in 1:length(A)], A.period; nperiod = A.nperiod)
    end
end
function LinearAlgebra.norm(A::PeriodicMatrix, p::Real = 2)
    n = norm([norm(A.M[i]) for i in 1:length(A)],p)
    if p == 2 
       return n*sqrt(A.nperiod)
    elseif isinf(p) 
       return n
    elseif p == 1 
       return n*A.nperiod
    else 
        throw(ArgumentError("only p-norms for p = 1, 2, or Inf are supported"))
    end
end
function +(A::PeriodicMatrix, B::PeriodicMatrix)
    isconstant(A) || isconstant(B) || A.Ts ≈ B.Ts || error("A and B must have the same sampling time")
    period = promote_period(A, B)
    pa = length(A) 
    pb = length(B)
    p = lcm(pa,pb)
    nta = numerator(rationalize(period/A.period))
    K = nta*A.nperiod*pa
    T = promote_type(eltype(A),eltype(B))
    X = Vector{Matrix{T}}(undef, p)
    for i = 1:p
        ia = mod(i-1,pa)+1
        ib = mod(i-1,pb)+1
        size(A.M[ia]) == size(B.M[ib]) || throw(DimensionMismatch("A and B have incompatible dimension"))
        X[i] = A.M[ia]+B.M[ib]
    end
    return PeriodicMatrix(X, period; nperiod = div(K,p))
end
+(A::PeriodicMatrix, C::AbstractMatrix) = +(A, PeriodicMatrix(C, A.Ts; nperiod = 1))
+(A::AbstractMatrix, C::PeriodicMatrix) = +(PeriodicMatrix(A, C.Ts; nperiod = 1), C)
-(A::PeriodicMatrix) = PeriodicMatrix(-A.M, A.period; nperiod = A.nperiod)
-(A::PeriodicMatrix, B::PeriodicMatrix) = +(A,-B)
-(A::PeriodicMatrix, C::AbstractMatrix) = +(A,-C)
-(A::AbstractMatrix, C::PeriodicMatrix) = +(A, -C)
function (+)(A::PeriodicMatrix, J::UniformScaling{<:Real})
    m, n = size(A)
    n == m || throw(DimensionMismatch("periodic matrix is not square: dimensions are $((m,n))"))
    PeriodicMatrix([A.M[i] + J for i in 1:length(A.M)], A.period; nperiod = A.nperiod)
end
(+)(J::UniformScaling{<:Real}, A::PeriodicMatrix) = +(A,J)
(-)(A::PeriodicMatrix, J::UniformScaling{<:Real}) = +(A,-J)
(-)(J::UniformScaling{<:Real}, A::PeriodicMatrix) = +(-A,J)
function pmsymadd!(A::Union{PeriodicMatrix,SwitchingPeriodicMatrix}, α = 1)
    # compute the symmetric matrix α*(A+transpose(A))
    m, n = size(A) 
    m == n || throw(ArgumentError("matrix A must be square"))
    if α == 1
       for i = 1:length(A) 
           inplace_transpose_add!(view(A.M[i],:,:))
       end
    else
        for i = 1:length(A) 
            inplace_transpose_add!(view(A.M[i],:,:),α)
        end
    end
    return A
end
function inplace_transpose_add!(A, α = 1)
    # compute (A+transpose(A))*α
    n = size(A, 1)
    if α == 1
        for i in 1:n
            for j in i:n
                A[i, j] += A[j, i]
                A[j, i] = A[i, j]
            end
        end
    else
        for i in 1:n
            for j in i:n
                A[i, j] += A[j, i]
                A[i,j]  *= α
                A[j, i] = A[i, j]
            end
        end
    end
end

function *(A::PeriodicMatrix, B::PeriodicMatrix)
    isconstant(A) || isconstant(B) || A.Ts ≈ B.Ts || error("A and B must have the same sampling time")
    period = promote_period(A, B)
    pa = length(A) 
    pb = length(B)
    nta = numerator(rationalize(period/A.period))
    ntb = numerator(rationalize(period/B.period))
    K = max(nta*A.nperiod*pa,ntb*B.nperiod*pb)
    p = lcm(pa,pb)
    T = promote_type(eltype(A),eltype(B))
    X = Vector{Matrix{T}}(undef, p)
    for i = 1:p
        ia = mod(i-1,pa)+1
        ib = mod(i-1,pb)+1
        size(A.M[ia],2) == size(B.M[ib],1) || throw(DimensionMismatch("A and B have incompatible dimension"))
        X[i] = A.M[ia]*B.M[ib]
    end
    return PeriodicMatrix(X, period; nperiod = div(K,p))
end
*(A::PeriodicMatrix, C::AbstractMatrix) = *(A, PeriodicMatrix(C, A.Ts; nperiod = 1))
*(A::AbstractMatrix, C::PeriodicMatrix) = *(PeriodicMatrix(A, C.Ts; nperiod = 1), C)
*(A::PeriodicMatrix, C::Real) = PeriodicMatrix(C.*A.M, A.period; nperiod = A.nperiod)
*(A::Real, C::PeriodicMatrix) = PeriodicMatrix(A.*C.M, C.period; nperiod = C.nperiod)
/(A::PeriodicMatrix, C::Real) = *(A, 1/C)
*(J::UniformScaling{<:Real}, A::PeriodicMatrix) = J.λ*A
*(A::PeriodicMatrix, J::UniformScaling{<:Real}) = A*J.λ
for (PMF, MF) in ((:pmmuladdsym, :muladdsym!), (:pmmultraddsym, :multraddsym!), (:pmmuladdtrsym,:muladdtrsym!) )
    @eval begin
        function $PMF(A::PeriodicMatrix,B::PeriodicMatrix,C::PeriodicMatrix, (α,β) = (true, true))
            # compute the symmetric results αA + βB*C, αA + βB'*C, and αA + βB*C'. 
            if isconstant(A)
               isconstant(B) || isconstant(C) || B.Ts ≈ C.Ts || error("B and C must have the same sampling time")
            elseif isconstant(B)
               isconstant(C) || A.Ts ≈ C.Ts || error("A and C must have the same sampling time")
            elseif isconstant(C)
                A.Ts ≈ B.Ts || error("A and B must have the same sampling time") 
            else
                A.Ts ≈ B.Ts ≈ C.Ts || error("A, B and C must have the same sampling time")
            end
            period = promote_period(A, B, C)
            pa = length(A) 
            pb = length(B)
            pc = length(C)
            p = lcm(pa,pb,pc)
            T = promote_type(eltype(A),eltype(B),eltype(B))
            X = Vector{Matrix{T}}(undef, p)
            nta = numerator(rationalize(period/A.period))
            ntb = numerator(rationalize(period/C.period))
            ntc = numerator(rationalize(period/C.period))
            Ka = nta*A.nperiod*pa
            Kb = ntb*B.nperiod*pb
            Kc = ntc*C.nperiod*pc
            K = max(Ka,Kb,Kc)
        
            for i = 1:p
                ia = mod(i-1,pa)+1
                ib = mod(i-1,pb)+1
                ic = mod(i-1,pc)+1
                # size(A.M[ia],2) == size(B.M[ib],1) || throw(DimensionMismatch("A and B have incompatible dimensions"))
                X[i] = copy(A.M[ia])
                $MF(view(X[i],:,:), B.M[ib], C.M[ic],(α,β))
            end
            return PeriodicMatrix(X, period; nperiod = div(K,p))
        end
        $PMF(A::PeriodicMatrix,B::PeriodicMatrix,C::PeriodicMatrix, α, β) = $PMF(A, B, C, (α,β))
    end
end
for (PMF, MF) in ((:pmmulsym, :muladdsym!), (:pmtrmulsym, :multraddsym!), (:pmmultrsym,:muladdtrsym!) )
    @eval begin
        function $PMF(B::PeriodicMatrix,C::PeriodicMatrix, β = true)
            # compute the symmetric results βB*C, βB'*C, and βB*C'. 
            isconstant(B) || isconstant(C) || B.Ts ≈ C.Ts || error("B and C must have the same sampling time")
            period = promote_period(B, C)
            pb = length(B)
            pc = length(C)
            p = lcm(pb,pc)
            T = promote_type(eltype(B),eltype(B))
            X = Vector{Matrix{T}}(undef, p)
            ntb = numerator(rationalize(period/B.period))
            ntc = numerator(rationalize(period/C.period))
            K = max(ntb*B.nperiod*pb,ntc*C.nperiod*pc)   
            tb = $PMF == pmtrmulsym
            for i = 1:p
                ib = mod(i-1,pb)+1
                ic = mod(i-1,pc)+1
                n = tb ? size(B.M[ib],2) : size(B.M[ib],1)
                X[i] = Matrix{T}(undef, n, n)
                $MF(view(X[i],:,:), B.M[ib], C.M[ic],(0,β))
            end
            return PeriodicMatrix(X, period; nperiod = div(K,p))
        end
    end
end

# muladdsym!(A::AbstractMatrix, B::AbstractMatrix, C::AbstractMatrix, α, β) = muladdsym!(A,B,C,(α,β))
function muladdsym!(A::AbstractMatrix, B::AbstractMatrix, C::AbstractMatrix, (α,β) = (true, true))
    # compute in A the symmetrix matrix α*A +  β*B*C
    n = LinearAlgebra.checksquare(A)
    n == size(B,1) || throw(ArgumentError("A and B must have the same number of rows"))
    n == size(C,2) || throw(ArgumentError("A and C must have the same number of columns"))
    m = size(B,2)
    m == size(C,1) || throw(ArgumentError("B and C have incompatible dimensions"))
    ZERO = zero(promote_type(eltype(B),eltype(C)))
    if α == 0
        for i = 1:n
            for j = i:n
                temp = ZERO
                for k = 1:m
                    temp += (B[i,k]*C[k,j])
                end
                A[i,j] = β*temp
                A[j,i] = A[i,j]
            end
        end
    else
        for i = 1:n
            for j = i:n
                temp = ZERO
                for k = 1:m
                    temp += (B[i,k]*C[k,j])
                end
                A[i,j] = α*A[i,j]+β*temp
                A[j,i] = A[i,j]
            end
        end
    end
    return A
end
# multraddsym!(A::AbstractMatrix, B::AbstractMatrix, C::AbstractMatrix, α, β) = multraddsym!(A, B, C, (α,β))
function multraddsym!(A::AbstractMatrix, B::AbstractMatrix, C::AbstractMatrix, (α,β) = (true, true))
    # compute in A the symmetrix matrix α*A +  β*transpose(B)*C
    n = LinearAlgebra.checksquare(A)
    n == size(B,2) || throw(ArgumentError("A and B must have the same number of columns"))
    n == size(C,2) || throw(ArgumentError("A and C must have the same number of columns"))
    m = size(B,1)
    m == size(C,1) || throw(ArgumentError("B and C have incompatible dimensions"))
    ZERO = zero(promote_type(eltype(B),eltype(C)))
    if α == 0
        for i = 1:n
            for j = i:n
                temp = ZERO
                for k = 1:m
                    temp += (B[k,i]*C[k,j])
                end
                A[i,j] = β*temp
                A[j,i] = A[i,j]
            end
        end
    else
        for i = 1:n
            for j = i:n
                temp = ZERO
                for k = 1:m
                    temp += (B[k,i]*C[k,j])
                end
                A[i,j] = α*A[i,j]+β*temp
                A[j,i] = A[i,j]
            end
        end
    end
    return A
end
# muladdtrsym!(A::AbstractMatrix, B::AbstractMatrix, C::AbstractMatrix, α, β) = muladdtrsym!(A, B, C, (α,β))
function muladdtrsym!(A::AbstractMatrix, B::AbstractMatrix, C::AbstractMatrix, (α,β) = (true, true))
    # compute in A the symmetrix matrix α*A +  β*transpose(B)*C
    n = LinearAlgebra.checksquare(A)
    n == size(B,1) || throw(ArgumentError("A and B must have the same number of rows"))
    n == size(C,1) || throw(ArgumentError("A and C must have the same number of rows"))
    m = size(B,2)
    m == size(C,2) || throw(ArgumentError("B and C have incompatible dimensions"))
    ZERO = zero(promote_type(eltype(B),eltype(C)))
    if α == 0
        for i = 1:n
            for j = i:n
                temp = ZERO
                for k = 1:m
                    temp += (B[i,k]*C[j,k])
                end
                A[i,j] = β*temp
                A[j,i] = A[i,j]
            end
        end
    else
        for i = 1:n
            for j = i:n
                temp = ZERO
                for k = 1:m
                    temp += (B[i,k]*C[j,k])
                end
                A[i,j] = α*A[i,j]+β*temp
                A[j,i] = A[i,j]
            end
        end
    end
    return A
end
for (PMF, MF) in ((:pmata, :mulatasym), (:pmaat, :mulaatsym) )
    @eval begin
        function $PMF(A::PeriodicMatrix)
            p = length(A) 
            T = eltype(A)
            X = Vector{Matrix{T}}(undef, p)
            for i = 1:p
                X[i] = $MF(view(A.M[i],:,:))
            end
            return PeriodicMatrix(X, A.period; nperiod = A.nperiod)
        end
    end
end
for (PMF, MF) in ((:pmata, :mulatasym), (:pmaat, :mulaatsym) )
    @eval begin
        function $PMF(A::PeriodicArray)
            m, n, p = size(A.M) 
            T = eltype(A)
            mn = $PMF == pmata ? n : m
            X = PeriodicArray(Array{T,3}(undef, mn, mn, p), A.period; nperiod = A.nperiod) 
            for i = 1:p
                copyto!(view(X.M,:,:,i), $MF(view(A.M,:,:,i)))
            end
            return X
        end
    end
end

"""
    pmata(A)

Compute for a discrete-time periodic matrix `A` the symmetric matrix `transpose(A)*A`. 

_Note:_ This function is available only for periodic matrices of types `PeriodicArray` and `PeriodicMatrix`.
"""
pmata(A)
"""
    pmaat(A)

Compute for a discrete-time periodic matrix `A` the symmetric matrix `A*transpose(A)`. 

_Note:_ This function is available only for periodic matrices of types `PeriodicArray` and `PeriodicMatrix`.
"""
pmaat(A)

        
function mulatasym(A::AbstractMatrix{T}) where {T}
    # compute the symmetrix matrix transpose(A)*A
    m, n = size(A,1), size(A,2) 
    X = similar(A,n,n)    
    ZERO = zero(T)
    for i = 1:n
        for j = i:n
            temp = ZERO
            for k = 1:m
                temp += (A[k,i]*A[k,j])
            end
            X[i,j] = temp
            X[j,i] = temp
        end
    end
    return X
end
function mulaatsym(A::AbstractMatrix{T}) where {T}
    # compute the symmetrix matrix A*transpose(A)
    m, n = size(A,1), size(A,2) 
    X = similar(A,m,m)    
    ZERO = zero(T)
    for i = 1:m
        for j = i:m
            temp = ZERO
            for k = 1:n
                temp += (A[i,k]*A[j,k])
            end
            X[i,j] = temp
            X[j,i] = temp
        end
    end
    return X
end


LinearAlgebra.issymmetric(A::PeriodicMatrix) = all([issymmetric(A.M[i]) for i in 1:A.dperiod])
Base.iszero(A::PeriodicMatrix) = iszero(A.M)
function ==(A::PeriodicMatrix, B::PeriodicMatrix)
    A.period == B.period && A.nperiod == B.nperiod && isequal(A.M, B.M) 
end
# function Base.isapprox(A::PeriodicMatrix, B::PeriodicMatrix; rtol::Real = sqrt(eps(Float64)), atol::Real = 0)
#     isconstant(A) && isconstant(B) && (return isapprox(A.M, B.M; rtol, atol))
#     isapprox(A.M, B.M; rtol, atol) && (A.period == B.period || A.period*B.nperiod == B.period*A.nperiod)
# end
function Base.isapprox(A::PeriodicMatrix, B::PeriodicMatrix; rtol::Real = sqrt(eps(Float64)), atol::Real = 0)
    A.period == B.period && A.nperiod == B.nperiod && isapprox(A.M, B.M; rtol, atol)
end
function Base.isapprox(A::PeriodicMatrix, J::UniformScaling{<:Real}; kwargs...)
    all([isapprox(A.M[i], J; kwargs...) for i in 1:length(A.M)])
end
Base.isapprox(J::UniformScaling{<:Real}, A::PeriodicMatrix; kwargs...) = isapprox(A, J; kwargs...)

function horzcat(A::PeriodicMatrix, B::PeriodicMatrix)
    A.Ts ≈ B.Ts || error("A and B must have the same sampling time")
    period = promote_period(A, B)
    pa = length(A) 
    pb = length(B)
    p = lcm(pa,pb)
    nta = numerator(rationalize(period/A.period))
    ntb = numerator(rationalize(period/B.period))
    K = max(nta*A.nperiod*pa,ntb*B.nperiod*pb)
    T = promote_type(eltype(A),eltype(B))
    X = Vector{Matrix{T}}(undef, p)
    for i = 1:p
        ia = mod(i-1,pa)+1
        ib = mod(i-1,pb)+1
        size(A.M[ia],1) == size(B.M[ib],1) || throw(DimensionMismatch("A and B have incompatible dimension"))
        X[i] = [A.M[ia] B.M[ib]]
    end
    return PeriodicMatrix(X, period; nperiod = div(K,p))
end
Base.hcat(A::PeriodicMatrix, B::PeriodicMatrix) = horzcat(A,B)
horzcat(A::PeriodicMatrix, B::AbstractMatrix) = horzcat(A, PeriodicMatrix(B, A.Ts; nperiod = 1))
horzcat(A::AbstractMatrix, B::PeriodicMatrix) = horzcat(PeriodicMatrix(A, B.Ts; nperiod = 1), B)
hcat(A::PeriodicMatrix, B::AbstractMatrix) = horzcat(A, PeriodicMatrix(B, A.Ts; nperiod = 1))
hcat(A::AbstractMatrix, B::PeriodicMatrix) = horzcat(PeriodicMatrix(A, B.Ts; nperiod = 1), B)

function vertcat(A::PeriodicMatrix, B::PeriodicMatrix)
    A.Ts ≈ B.Ts || error("A and B must have the same sampling time")
    period = promote_period(A, B)
    pa = length(A) 
    pb = length(B)
    p = lcm(pa,pb)
    nta = numerator(rationalize(period/A.period))
    ntb = numerator(rationalize(period/B.period))
    K = max(nta*A.nperiod*pa,ntb*B.nperiod*pb)
    T = promote_type(eltype(A),eltype(B))
    X = Vector{Matrix{T}}(undef, p)
    for i = 1:p
        ia = mod(i-1,pa)+1
        ib = mod(i-1,pb)+1
        size(A.M[ia],2) == size(B.M[ib],2) || throw(DimensionMismatch("A and B have incompatible dimension"))
        X[i] = [A.M[ia]; B.M[ib]]
    end
    return PeriodicMatrix(X, period; nperiod = div(K,p))
end
Base.vcat(A::PeriodicMatrix, B::PeriodicMatrix) = vertcat(A,B)
vertcat(A::PeriodicMatrix, B::AbstractMatrix) = vertcat(A, PeriodicMatrix(B, A.Ts; nperiod = 1))
vertcat(A::AbstractMatrix, B::PeriodicMatrix) = vertcat(PeriodicMatrix(A, B.Ts; nperiod = 1), B)
vcat(A::PeriodicMatrix, B::AbstractMatrix) = vertcat(A, PeriodicMatrix(B, A.Ts; nperiod = 1))
vcat(A::AbstractMatrix, B::PeriodicMatrix) = vertcat(PeriodicMatrix(A, B.Ts; nperiod = 1), B)

function blockdiag(A::PeriodicMatrix, B::PeriodicMatrix)
    A.Ts ≈ B.Ts || error("A and B must have the same sampling time")
    period = promote_period(A, B)
    pa = length(A) 
    pb = length(B)
    p = lcm(pa,pb)
    nta = numerator(rationalize(period/A.period))
    ntb = numerator(rationalize(period/B.period))
    K = max(nta*A.nperiod*pa,ntb*B.nperiod*pb)
    T = promote_type(eltype(A),eltype(B))
    X = Vector{Matrix{T}}(undef, p)
    for i = 1:p
        ia = mod(i-1,pa)+1
        ib = mod(i-1,pb)+1
        X[i] = bldiag(A.M[ia], B.M[ib])
    end
    return PeriodicMatrix(X, period; nperiod = div(K,p))
end


#  SwitchingPeriodicMatrix
function pmshift(A::SwitchingPeriodicMatrix, k::Int = 1)
    return convert(SwitchingPeriodicMatrix,pmshift(convert(PeriodicMatrix,A),k))
end
LinearAlgebra.inv(A::SwitchingPeriodicMatrix) = SwitchingPeriodicMatrix(inv.(A.M), A.ns, A.period; nperiod = A.nperiod)
function LinearAlgebra.transpose(A::SwitchingPeriodicMatrix)
    return SwitchingPeriodicMatrix(copy.(transpose.(A.M)), A.ns, A.period; nperiod = A.nperiod)
end
function LinearAlgebra.adjoint(A::SwitchingPeriodicMatrix)
    return SwitchingPeriodicMatrix(copy.(transpose.(A.M)), A.ns, A.period; nperiod = A.nperiod)
end
function Base.reverse(A::SwitchingPeriodicMatrix)
    n = length(A.ns)
    return SwitchingPeriodicMatrix(reverse(A.M), n == 1 ? A.ns : [A.ns[n].-reverse(A.ns[1:n-1]); A.ns[n]], A.period; nperiod = A.nperiod)
end
function LinearAlgebra.opnorm(A::SwitchingPeriodicMatrix, p::Union{Real,Missing} = missing)
    if ismissing(p)
       return SwitchingPeriodicMatrix([norm(A.M[i]) for i in 1:length(A.M)], A.ns, A.period; nperiod = A.nperiod)
    else
       return SwitchingPeriodicMatrix([opnorm(A.M[i],p) for i in 1:length(A.M)], A.ns, A.period; nperiod = A.nperiod)
    end
end
function LinearAlgebra.norm(A::SwitchingPeriodicMatrix, p::Real = 2)
    k = length(A)
    k == 0 && (return zero(eltype(A)))
    tn = norm.(A.M)
    if p == 2
       tn[1] *= sqrt(A.ns[1]) 
       for i = 2:k
           tn[i] *= sqrt(A.ns[i]-A.ns[i-1])
       end
       return norm(tn,p)*sqrt(A.nperiod)
    elseif p == 1
       tn[1] *= A.ns[1] 
       for i = 2:k
           tn[i] *= (A.ns[i]-A.ns[i-1])
       end
       return norm(tn,p)*A.nperiod
    else 
       isinf(p) || throw(ArgumentError("only p-norms for p = 1, 2, or Inf are supported"))
       return norm(tn,p) 
    end
end
function LinearAlgebra.tr(A::SwitchingPeriodicMatrix)
    return SwitchingPeriodicMatrix([[tr(A.M[i]);;] for i in 1:length(A.M)], A.ns, A.period; nperiod = A.nperiod)
end
function trace(A::SwitchingPeriodicMatrix) 
    t = zero(eltype(A))
    k = length(A)
    k == 0 && (return t)
    t += tr(A.M[1])*A.ns[1]
    for i in 2:k
        t += tr(A.M[i])*(A.ns[i]-A.ns[i-1])
    end
    return t*A.nperiod
end

Base.iszero(A::SwitchingPeriodicMatrix) = iszero(A.M)
LinearAlgebra.issymmetric(A::SwitchingPeriodicMatrix) = all([issymmetric(A.M[i]) for i in 1:length(A.ns)])
function ==(A::SwitchingPeriodicMatrix, B::SwitchingPeriodicMatrix)
    A.period == B.period && A.nperiod == B.nperiod && A.ns == B.ns && isequal(A.M, B.M) 
end
function Base.isapprox(A::SwitchingPeriodicMatrix, B::SwitchingPeriodicMatrix; rtol::Real = sqrt(eps(Float64)), atol::Real = 0)
    A.period == B.period && A.nperiod == B.nperiod && A.ns == B.ns && isapprox(A.M, B.M; rtol, atol)
end
function Base.isapprox(A::SwitchingPeriodicMatrix, J::UniformScaling{<:Real}; kwargs...)
    all([isapprox(A.M[i], J; kwargs...) for i in 1:length(A.M)])
end
Base.isapprox(J::UniformScaling{<:Real}, A::SwitchingPeriodicMatrix; kwargs...) = isapprox(A, J; kwargs...)

function +(A::SwitchingPeriodicMatrix, B::SwitchingPeriodicMatrix)
    A.Ts ≈ B.Ts || error("A and B must have the same sampling time")
    A.period == B.period && A.nperiod == B.nperiod && A.ns == B.ns &&
    (return SwitchingPeriodicMatrix{:d,promote_type(eltype(A),eltype(B))}([A.M[i]+B.M[i] for i in 1:length(A.M)], A.ns, A.period, A.nperiod))
    isconstant(A) && 
       (return SwitchingPeriodicMatrix{:d,promote_type(eltype(A),eltype(B))}([A.M[1]+B.M[i] for i in 1:length(B.M)], B.ns, B.period, B.nperiod))
    isconstant(B) && 
       (return SwitchingPeriodicMatrix{:d,promote_type(eltype(A),eltype(B))}([A.M[i]+B.M[1] for i in 1:length(A.M)], A.ns, A.period, A.nperiod))
    if A.period == B.period 
       nperiod = A.nperiod
       if nperiod == B.nperiod
          ns = unique(sort([A.ns;B.ns]))
       else
          nperiod = gcd(A.nperiod,B.nperiod)
          ns = unique(sort([vcat([(i-1)*A.dperiod .+ A.ns for i in 1:div(A.nperiod,nperiod)]...);
                          vcat([(i-1)*B.dperiod .+ B.ns for i in 1:div(B.nperiod,nperiod)]...)]))
       end
       N = length(ns)                  
       return SwitchingPeriodicMatrix{:d,promote_type(eltype(A),eltype(B))}([kpmeval(A,ns[i])+kpmeval(B,ns[i]) for i in 1:N], ns, A.period, nperiod)   
    else
        # TODO: implement explicit computations
        period = promote_period(A,B)
        return set_period(A,period)+set_period(B,period)
        #error("periods must be equal for addition")
    end
end
+(A::SwitchingPeriodicMatrix, C::AbstractMatrix) = +(A, SwitchingPeriodicMatrix([C], [A.dperiod], A.period))
+(A::AbstractMatrix, C::SwitchingPeriodicMatrix) = +(SwitchingPeriodicMatrix([A], [C.dperiod], C.period), C)
-(A::SwitchingPeriodicMatrix) = SwitchingPeriodicMatrix{:d,eltype(A)}(-A.M, A.ns, A.period, A.nperiod)
-(A::SwitchingPeriodicMatrix, B::SwitchingPeriodicMatrix) = +(A,-B)
-(A::SwitchingPeriodicMatrix, C::AbstractMatrix) = +(A,-C)
-(A::AbstractMatrix, C::SwitchingPeriodicMatrix) = +(A, -C)
function (+)(A::SwitchingPeriodicMatrix, J::UniformScaling{<:Real}) 
    mv, nv = size(A)
    n = minimum(nv);  
    n == maximum(nv) == minimum(mv) == maximum(mv) || throw(DimensionMismatch("periodic matrix is not square: dimensions are $((mv,nv))"))
    A+Matrix(J(n))
end
(+)(J::UniformScaling{<:Real}, A::SwitchingPeriodicMatrix) = +(A,J)
(-)(A::SwitchingPeriodicMatrix, J::UniformScaling{<:Real}) = +(A,-J)
(-)(J::UniformScaling{<:Real}, A::SwitchingPeriodicMatrix) = +(-A,J)

function *(A::SwitchingPeriodicMatrix, B::SwitchingPeriodicMatrix)
    A.Ts ≈ B.Ts || error("A and B must have the same sampling time")
    A.period == B.period && A.nperiod == B.nperiod && A.ns == B.ns &&
    (return SwitchingPeriodicMatrix{:d,promote_type(eltype(A),eltype(B))}([A.M[i]*B.M[i] for i in 1:length(A.M)], A.ns, A.period, A.nperiod))
    isconstant(A) && 
       (return SwitchingPeriodicMatrix{:d,promote_type(eltype(A),eltype(B))}([A.M[1]*B.M[i] for i in 1:length(B.M)], B.ns, B.period, B.nperiod))
    isconstant(B) && 
       (return SwitchingPeriodicMatrix{:d,promote_type(eltype(A),eltype(B))}([A.M[i]*B.M[1] for i in 1:length(A.M)], A.ns, A.period, A.nperiod))
    if A.period == B.period 
       nperiod = A.nperiod
       if nperiod == B.nperiod
          ns = unique(sort([A.ns;B.ns]))
       else
          nperiod = gcd(A.nperiod,B.nperiod)
          ns = unique(sort([vcat([(i-1)*A.dperiod .+ A.ns for i in 1:div(A.nperiod,nperiod)]...);
                          vcat([(i-1)*B.dperiod .+ B.ns for i in 1:div(B.nperiod,nperiod)]...)]))
       end
       N = length(ns)                  
       return SwitchingPeriodicMatrix{:d,promote_type(eltype(A),eltype(B))}([kpmeval(A,ns[i])*kpmeval(B,ns[i]) for i in 1:N], ns, A.period, nperiod)   
    else
       # TODO: implement explicit computations
       period = promote_period(A,B)
       return set_period(A,period)*set_period(B,period)
       #error("periods must be equal for multiplication")
    end
end
*(A::SwitchingPeriodicMatrix, C::AbstractMatrix) = SwitchingPeriodicMatrix{:d,promote_type(eltype(A),eltype(C))}([A.M[i]*C for i in 1:length(A.M)], A.ns, A.period, A.nperiod)
*(A::AbstractMatrix, C::SwitchingPeriodicMatrix) = SwitchingPeriodicMatrix{:d,promote_type(eltype(A),eltype(C))}([A*C.M[i] for i in 1:length(C.M)], C.ns, C.period, C.nperiod)
*(A::SwitchingPeriodicMatrix, C::Real) = SwitchingPeriodicMatrix{:d,promote_type(eltype(A),eltype(C))}([A.M[i]*C for i in 1:length(A.M)], A.ns, A.period, A.nperiod)
*(C::Real, A::SwitchingPeriodicMatrix) = SwitchingPeriodicMatrix{:d,promote_type(eltype(A),eltype(C))}([C*A.M[i] for i in 1:length(A.M)], A.ns, A.period, A.nperiod)
/(A::SwitchingPeriodicMatrix, C::Real) = *(A, 1/C)
*(J::UniformScaling{<:Real}, A::SwitchingPeriodicMatrix) = J.λ*A
*(A::SwitchingPeriodicMatrix, J::UniformScaling{<:Real}) = A*J.λ

function horzcat(A::SwitchingPeriodicMatrix, B::SwitchingPeriodicMatrix)
    A.period == B.period && A.nperiod == B.nperiod && A.ns == B.ns &&
        (return SwitchingPeriodicMatrix{:d,promote_type(eltype(A),eltype(B))}([[A.M[i] B.M[i]] for i in 1:length(A.M)], A.ns, A.period, A.nperiod))
    isconstant(A) && 
       (return SwitchingPeriodicMatrix{:d,promote_type(eltype(A),eltype(B))}([[A.M[1] B.M[i]] for i in 1:length(B.M)], B.ns, B.period, B.nperiod))
    isconstant(B) && 
       (return SwitchingPeriodicMatrix{:d,promote_type(eltype(A),eltype(B))}([[A.M[i] B.M[1]] for i in 1:length(A.M)], A.ns, A.period, A.nperiod))
    if A.period == B.period 
       nperiod = A.nperiod
       if  nperiod == B.nperiod
           ns = unique(sort([A.ns;B.ns]))
       else
           nperiod = gcd(A.nperiod,B.nperiod)
           ns = unique(sort([vcat([(i-1)*A.dperiod .+ A.ns for i in 1:div(A.nperiod,nperiod)]...);
                          vcat([(i-1)*B.dperiod .+ B.ns for i in 1:div(B.nperiod,nperiod)]...)]))
       end
       N = length(ns)                  
       return SwitchingPeriodicMatrix{:d,promote_type(eltype(A),eltype(B))}([[kpmeval(A,ns[i]) kpmeval(B,ns[i])] for i in 1:N], ns, A.period, nperiod)   
    else 
        # TODO: implement explicit computations
        period = promote_period(A,B)
        return horzcat(set_period(A,period),set_period(B,period))
        #error("periods must be equal for horizontal concatenation")
    end
end
Base.hcat(A::SwitchingPeriodicMatrix, B::SwitchingPeriodicMatrix) = horzcat(A,B)
horzcat(A::SwitchingPeriodicMatrix, B::AbstractMatrix) = horzcat(A, SwitchingPeriodicMatrix([B], [A.dperiod], A.period))
horzcat(A::AbstractMatrix, B::SwitchingPeriodicMatrix) = horzcat(SwitchingPeriodicMatrix([A], [B.dperiod], B.period), B)
hcat(A::SwitchingPeriodicMatrix, B::AbstractMatrix) = horzcat(A, SwitchingPeriodicMatrix([B], [A.dperiod], A.period))
hcat(A::AbstractMatrix, B::SwitchingPeriodicMatrix) = horzcat(SwitchingPeriodicMatrix([A], [B.dperiod], B.period), B)


function vertcat(A::SwitchingPeriodicMatrix, B::SwitchingPeriodicMatrix)
    A.period == B.period && A.nperiod == B.nperiod && A.ns == B.ns &&
        (return SwitchingPeriodicMatrix{:d,promote_type(eltype(A),eltype(B))}([[A.M[i]; B.M[i]] for i in 1:length(A.M)], A.ns, A.period, A.nperiod))
    isconstant(A) && 
       (return SwitchingPeriodicMatrix{:d,promote_type(eltype(A),eltype(B))}([[A.M[1]; B.M[i]] for i in 1:length(B.M)], B.ns, B.period, B.nperiod))
    isconstant(B) && 
       (return SwitchingPeriodicMatrix{:d,promote_type(eltype(A),eltype(B))}([[A.M[i]; B.M[1]] for i in 1:length(A.M)], A.ns, A.period, A.nperiod))
    if A.period == B.period 
       nperiod = A.nperiod
       if  nperiod == B.nperiod
           ns = unique(sort([A.ns;B.ns]))
       else
           nperiod = gcd(A.nperiod,B.nperiod)
           ns = unique(sort([vcat([(i-1)*A.dperiod .+ A.ns for i in 1:div(A.nperiod,nperiod)]...);
                          vcat([(i-1)*B.dperiod .+ B.ns for i in 1:div(B.nperiod,nperiod)]...)]))
       end
       N = length(ns)                  
       return SwitchingPeriodicMatrix{:d,promote_type(eltype(A),eltype(B))}([[kpmeval(A,ns[i]); kpmeval(B,ns[i])] for i in 1:N], ns, A.period, nperiod)   
    else
        # TODO: implement explicit computations
        period = promote_period(A,B)
        return vertcat(set_period(A,period),set_period(B,period))
        #error("periods must be equal for vertical concatenation")
    end
end
Base.vcat(A::SwitchingPeriodicMatrix, B::SwitchingPeriodicMatrix) = vertcat(A,B)
vertcat(A::SwitchingPeriodicMatrix, B::AbstractMatrix) = vertcat(A, SwitchingPeriodicMatrix([B], [A.dperiod], A.period))
vertcat(A::AbstractMatrix, B::SwitchingPeriodicMatrix) = vertcat(SwitchingPeriodicMatrix([A], [B.dperiod], B.period), B)
vcat(A::SwitchingPeriodicMatrix, B::AbstractMatrix) = vertcat(A, SwitchingPeriodicMatrix([B], [A.dperiod], A.period))
vcat(A::AbstractMatrix, B::SwitchingPeriodicMatrix) = vertcat(SwitchingPeriodicMatrix([A], [B.dperiod], B.period), B)

function blockdiag(A::SwitchingPeriodicMatrix, B::SwitchingPeriodicMatrix)
    A.period == B.period && A.nperiod == B.nperiod && A.ns == B.ns &&
        (return SwitchingPeriodicMatrix{:d,promote_type(eltype(A),eltype(B))}([bldiag(A.M[i], B.M[i]) for i in 1:length(A.M)], A.ns, A.period, A.nperiod))
    isconstant(A) && 
       (return SwitchingPeriodicMatrix{:d,promote_type(eltype(A),eltype(B))}([bldiag(A.M[1], B.M[i]) for i in 1:length(B.M)], B.ns, B.period, B.nperiod))
    isconstant(B) && 
       (return SwitchingPeriodicMatrix{:d,promote_type(eltype(A),eltype(B))}([bldiag(A.M[i], B.M[1]) for i in 1:length(A.M)], A.ns, A.period, A.nperiod))

    if A.period == B.period 
       nperiod = A.nperiod
       if  nperiod == B.nperiod
           ns = unique(sort([A.ns;B.ns]))
       else
           nperiod = gcd(A.nperiod,B.nperiod)
           ns = unique(sort([vcat([(i-1)*A.dperiod .+ A.ns for i in 1:div(A.nperiod,nperiod)]...);
                          vcat([(i-1)*B.dperiod .+ B.ns for i in 1:div(B.nperiod,nperiod)]...)]))
       end
       N = length(ns)                  
       return SwitchingPeriodicMatrix{:d,promote_type(eltype(A),eltype(B))}([bldiag(kpmeval(A,ns[i]), kpmeval(B,ns[i])) for i in 1:N], ns, A.period, nperiod)   
    else
        # TODO: implement explicit computations
        period = promote_period(A,B)
        return blockdiag(set_period(A,period),set_period(B,period))
        #error("periods must be equal for bloc diagonal appending")
    end
end

#  SwitchingPeriodicArray
function pmshift(A::SwitchingPeriodicArray, k::Int = 1)
    return convert(SwitchingPeriodicArray,pmshift(convert(PeriodicArray,A),k))
end
function LinearAlgebra.inv(A::SwitchingPeriodicArray)
    x = similar(A.M)
    [x[:,:,i] = inv(A.M[:,:,i]) for i in 1:size(A.M,3)]
    SwitchingPeriodicArray(x, A.ns, A.period; nperiod = A.nperiod)
end
function LinearAlgebra.transpose(A::SwitchingPeriodicArray)
    return SwitchingPeriodicArray(permutedims(A.M,(2,1,3)), A.ns, A.period; nperiod = A.nperiod)
end
function LinearAlgebra.adjoint(A::SwitchingPeriodicArray)
    return SwitchingPeriodicArray(permutedims(A.M,(2,1,3)), A.ns, A.period; nperiod = A.nperiod)
end
function Base.reverse(A::SwitchingPeriodicArray)
    n = length(A.ns)
    return SwitchingPeriodicArray(reverse(A.M,dims=3), n == 1 ? A.ns : [A.ns[n].-reverse(A.ns[1:n-1]); A.ns[n]], A.period; nperiod = A.nperiod)
end
function LinearAlgebra.opnorm(A::SwitchingPeriodicArray, p::Union{Real,Missing} = missing)
    k = size(A.M,3)
    x = Array{eltype(A),3}(undef, 1, 1, k)
    if ismissing(p)
        [x[1,1,i] = norm(view(A.M,:,:,i)) for i in 1:k]  # Frobenius noorm
    else
        [x[1,1,i] = opnorm(view(A.M,:,:,i),p) for i in 1:k] # p-norm
    end
    return SwitchingPeriodicArray(x, A.ns, A.period; nperiod = A.nperiod)
end
function LinearAlgebra.norm(A::SwitchingPeriodicArray, p::Real = 2)
    k = length(A)
    k == 0 && (return zero(eltype(A)))
    tn = [norm(view(A.M,:,:,i)) for i in 1:k]
    if p == 2
       tn[1] *= sqrt(A.ns[1]) 
       for i = 2:k
           tn[i] *= sqrt(A.ns[i]-A.ns[i-1])
       end
       return norm(tn,p)*sqrt(A.nperiod)
    elseif p == 1
       tn[1] *= A.ns[1] 
       for i = 2:k
           tn[i] *= (A.ns[i]-A.ns[i-1])
       end
       return norm(tn,p)*A.nperiod
    else 
       isinf(p) || throw(ArgumentError("only p-norms for p = 1, 2, or Inf are supported"))
       return norm(tn,p) 
    end
end
function LinearAlgebra.tr(A::SwitchingPeriodicArray)
    p = size(A.M,3)
    x = Array{eltype(A),3}(undef, 1, 1, p)
    [x[1,1,i] = tr(view(A.M,:,:,i)) for i in 1:p]
    return SwitchingPeriodicArray(x, A.ns, A.period; nperiod = A.nperiod)
end
function trace(A::SwitchingPeriodicArray) 
    t = zero(eltype(A))
    k = length(A)
    k == 0 && (return t)
    t += tr(view(A.M,:,:,1))*A.ns[1]
    for i in 2:k
        t += tr(view(A.M,:,:,i))*(A.ns[i]-A.ns[i-1])
    end
    return t*A.nperiod
end

LinearAlgebra.issymmetric(A::SwitchingPeriodicArray) = all([issymmetric(view(A.M,:,:,i)) for i in 1:length(A.ns)])
Base.iszero(A::SwitchingPeriodicArray) = iszero(A.M)
function ==(A::SwitchingPeriodicArray, B::SwitchingPeriodicArray)
    A.period == B.period && A.nperiod == B.nperiod && A.ns == B.ns && isequal(A.M, B.M) 
end
function Base.isapprox(A::SwitchingPeriodicArray, B::SwitchingPeriodicArray; rtol::Real = sqrt(eps(Float64)), atol::Real = 0)
    A.period == B.period && A.nperiod == B.nperiod && A.ns == B.ns &&  isapprox(A.M, B.M; rtol, atol)
end
function Base.isapprox(A::SwitchingPeriodicArray, J::UniformScaling{<:Real}; kwargs...)
    all([isapprox(view(A.M,:,:,i), J; kwargs...) for i in 1:length(A.ns)])
end
Base.isapprox(J::UniformScaling{<:Real}, A::SwitchingPeriodicArray; kwargs...) = isapprox(A, J; kwargs...)

function +(A::SwitchingPeriodicArray, B::SwitchingPeriodicArray)
    isconstant(A) || isconstant(B) || A.Ts ≈ B.Ts || error("A and B must have the same sampling time")
    period = promote_period(A, B)
    T = promote_type(eltype(A),eltype(B))
    m, n, pa = size(A.M)
    mb, nb, pb = size(B.M)
    (m, n) == (mb, nb) || throw(DimensionMismatch("A and B must have the same dimensions"))

    A.period == B.period && A.nperiod == B.nperiod && A.ns == B.ns &&
    (return SwitchingPeriodicArray{:d,T}(A.M+B.M, A.ns, A.period, A.nperiod))
    if isconstant(A) 
        N = length(B.ns)
        X = Array{T,3}(undef, m, n, N)
        for i = 1:N
            X[:,:,i] = A.M[:,:,1]+B.M[:,:,i]
        end
        return SwitchingPeriodicArray{:d,T}(X, B.ns, B.period, B.nperiod)
    elseif isconstant(B) 
        N = length(A.ns)
        X = Array{T,3}(undef, m, n, N)
        for i = 1:N
            X[:,:,i] = A.M[:,:,i]+B.M[:,:,1]
        end
        return SwitchingPeriodicArray{:d,T}(X, A.ns, A.period, A.nperiod)
    elseif A.period == B.period 
        nperiod = A.nperiod
        if nperiod == B.nperiod
           ns = unique(sort([A.ns;B.ns]))
        else
           nperiod = gcd(A.nperiod,B.nperiod)
           ns = unique(sort([vcat([(i-1)*A.dperiod .+ A.ns for i in 1:div(A.nperiod,nperiod)]...);
                             vcat([(i-1)*B.dperiod .+ B.ns for i in 1:div(B.nperiod,nperiod)]...)]))
        end
        N = length(ns)   
        X = Array{T,3}(undef, m, n, N)
        for i = 1:N
            X[:,:,i] = kpmeval(A,ns[i])+kpmeval(B,ns[i])
        end
        return SwitchingPeriodicArray{:d,promote_type(eltype(A),eltype(B))}(X, ns, A.period, nperiod)   
    else
        # TODO: implement explicit computations
        period = promote_period(A,B)
        return set_period(A,period)+set_period(B,period)
        #error("periods must be equal for addition")
    end              
end
+(A::SwitchingPeriodicArray, C::AbstractMatrix) = +(A, SwitchingPeriodicArray(reshape(C,size(C,1),size(C,2),1), [A.dperiod], A.period))
+(A::AbstractMatrix, C::SwitchingPeriodicArray) = +(SwitchingPeriodicArray(reshape(A,size(A,1),size(A,2),1), [C.dperiod], C.period), C)
-(A::SwitchingPeriodicArray) = SwitchingPeriodicArray{:d,eltype(A)}(-A.M, A.ns, A.period, A.nperiod)
-(A::SwitchingPeriodicArray, B::SwitchingPeriodicArray) = +(A,-B)
-(A::SwitchingPeriodicArray, C::AbstractMatrix) = +(A,-C)
-(A::AbstractMatrix, C::SwitchingPeriodicArray) = +(A, -C)
function (+)(A::SwitchingPeriodicArray, J::UniformScaling{<:Real}) 
    m, n = size(A)
    n == m || throw(DimensionMismatch("periodic matrix is not square: dimensions are $((m,n))"))
    A+Matrix(J(n))
end
(+)(J::UniformScaling{<:Real}, A::SwitchingPeriodicArray) = +(A,J)
(-)(A::SwitchingPeriodicArray, J::UniformScaling{<:Real}) = +(A,-J)
(-)(J::UniformScaling{<:Real}, A::SwitchingPeriodicArray) = +(-A,J)

function *(A::SwitchingPeriodicArray, B::SwitchingPeriodicArray)
    isconstant(A) || isconstant(B) || A.Ts ≈ B.Ts || error("A and B must have the same sampling time")
    period = promote_period(A, B)
    T = promote_type(eltype(A),eltype(B))
    m, n, pa = size(A.M)
    mb, nb, pb = size(B.M)
    n == mb || throw(DimensionMismatch("number of columns of A $n not equal to the number of rows of B $mb"))

    if A.period == B.period && A.nperiod == B.nperiod && A.ns == B.ns 
        N = length(A.ns)
        X = Array{T,3}(undef, m, nb, N)
        for i = 1:N
            mul!(view(X,:,:,i),view(A.M,:,:,i),view(B.M,:,:,i))
        end
        return SwitchingPeriodicArray{:d,T}(X, A.ns, A.period, A.nperiod)
    end
    if isconstant(A) 
        N = length(B.ns)
        X = Array{T,3}(undef, m, nb, N)
        for i = 1:N
            mul!(view(X,:,:,i),view(A.M,:,:,1),view(B.M,:,:,i))
        end
        return SwitchingPeriodicArray{:d,T}(X, B.ns, B.period, B.nperiod)
    elseif isconstant(B) 
        N = length(A.ns)
        X = Array{T,3}(undef, m, n, N)
        for i = 1:N
            mul!(view(X,:,:,i),view(A.M,:,:,i),view(B.M,:,:,1))
        end
        return SwitchingPeriodicArray{:d,T}(X, A.ns, A.period, A.nperiod)
    elseif A.period == B.period 
        nperiod = A.nperiod
        if nperiod == B.nperiod
           ns = unique(sort([A.ns;B.ns]))
        else
           nperiod = gcd(A.nperiod,B.nperiod)
           ns = unique(sort([vcat([(i-1)*A.dperiod .+ A.ns for i in 1:div(A.nperiod,nperiod)]...);
                             vcat([(i-1)*B.dperiod .+ B.ns for i in 1:div(B.nperiod,nperiod)]...)]))
        end
        N = length(ns)   
        X = Array{T,3}(undef, m, nb, N)
        for i = 1:N
            mul!(view(X,:,:,i),kpmeval(A,ns[i]),kpmeval(B,ns[i]))
        end
        return SwitchingPeriodicArray{:d,promote_type(eltype(A),eltype(B))}(X, ns, A.period, nperiod)   
    else
        # TODO: implement explicit computations
        period = promote_period(A,B)
        return set_period(A,period)*set_period(B,period)
        #error("periods must be equal for multiplication")
    end
end
*(A::SwitchingPeriodicArray, C::AbstractMatrix) = *(A, SwitchingPeriodicArray(reshape(C,size(C,1),size(C,2),1), [A.dperiod], A.period))
*(A::AbstractMatrix, C::SwitchingPeriodicArray) = *(SwitchingPeriodicArray(reshape(A,size(A,1),size(A,2),1), [C.dperiod], C.period), C)
# *(A::SwitchingPeriodicArray, C::AbstractMatrix) = SwitchingPeriodicArray{:d,promote_type(eltype(A),eltype(C))}([A.M[i]*C for i in 1:length(A.M)], A.ns, A.period, A.nperiod)
# *(A::AbstractMatrix, C::SwitchingPeriodicArray) = SwitchingPeriodicArray{:d,promote_type(eltype(A),eltype(C))}([A*C.M[i] for i in 1:length(C.M)], C.ns, C.period, C.nperiod)
*(A::SwitchingPeriodicArray, C::Real) = SwitchingPeriodicArray{:d,promote_type(eltype(A),eltype(C))}(A.M*C, A.ns, A.period, A.nperiod)
*(C::Real, A::SwitchingPeriodicArray) = SwitchingPeriodicArray{:d,promote_type(eltype(A),eltype(C))}(C*A.M, A.ns, A.period, A.nperiod)
/(A::SwitchingPeriodicArray, C::Real) = *(A, 1/C)
*(J::UniformScaling{<:Real}, A::SwitchingPeriodicArray) = J.λ*A
*(A::SwitchingPeriodicArray, J::UniformScaling{<:Real}) = A*J.λ


function horzcat(A::SwitchingPeriodicArray, B::SwitchingPeriodicArray)
    isconstant(A) || isconstant(B) || A.Ts ≈ B.Ts || error("A and B must have the same sampling time")
    period = promote_period(A, B)
    T = promote_type(eltype(A),eltype(B))
    ma, na, pa = size(A.M)
    mb, nb, pb = size(B.M)
    ma == mb || throw(DimensionMismatch("number of rows of A $ma not equal to the number of rows of B $mb"))

    if A.period == B.period && A.nperiod == B.nperiod && A.ns == B.ns 
        N = length(A.ns)
        X = Array{T,3}(undef, ma, na+nb, N)
        j1 = 1:na; j2 = na+1:na+nb
        for i = 1:N
            copyto!(view(X,:,j1,i),view(A.M,:,:,i))
            copyto!(view(X,:,j2,i),view(B.M,:,:,i))
        end
        return SwitchingPeriodicArray{:d,T}(X, A.ns, A.period, A.nperiod)
    end
    if isconstant(A) 
        N = length(B.ns)
        X = Array{T,3}(undef, ma, na+nb, N)
        j1 = 1:na; j2 = na+1:na+nb
        for i = 1:N
            copyto!(view(X,:,j1,i),view(A.M,:,:,1))
            copyto!(view(X,:,j2,i),view(B.M,:,:,i))
        end
        return SwitchingPeriodicArray{:d,T}(X, B.ns, B.period, B.nperiod)
    elseif isconstant(B) 
        N = length(A.ns)
        X = Array{T,3}(undef, ma, na+nb, N)
        j1 = 1:na; j2 = na+1:na+nb
        for i = 1:N
            copyto!(view(X,:,j1,i),view(A.M,:,:,i))
            copyto!(view(X,:,j2,i),view(B.M,:,:,1))
        end
        return SwitchingPeriodicArray{:d,T}(X, A.ns, A.period, A.nperiod)
    elseif A.period == B.period 
        nperiod = A.nperiod
        if nperiod == B.nperiod
           ns = unique(sort([A.ns;B.ns]))
        else
           nperiod = gcd(A.nperiod,B.nperiod)
           ns = unique(sort([vcat([(i-1)*A.dperiod .+ A.ns for i in 1:div(A.nperiod,nperiod)]...);
                             vcat([(i-1)*B.dperiod .+ B.ns for i in 1:div(B.nperiod,nperiod)]...)]))
        end
        N = length(ns)   
        X = Array{T,3}(undef, ma, na+nb, N)
        j1 = 1:na; j2 = na+1:na+nb
        for i = 1:N
            copyto!(view(X,:,j1,i),kpmeval(A,ns[i]))
            copyto!(view(X,:,j2,i),kpmeval(B,ns[i]))
        end
        return SwitchingPeriodicArray{:d,promote_type(eltype(A),eltype(B))}(X, ns, A.period, nperiod)  
    else 
        # TODO: implement explicit computations
        period = promote_period(A,B)
        return horzcat(set_period(A,period),set_period(B,period))
        #error("periods must be equal for horizontal concatenation")
    end
end
Base.hcat(A::SwitchingPeriodicArray, B::SwitchingPeriodicArray) = horzcat(A,B)
horzcat(A::SwitchingPeriodicArray, B::AbstractMatrix) = horzcat(A, SwitchingPeriodicArray(reshape(B,size(B,1),size(B,2),1), [A.dperiod], A.period))
horzcat(A::AbstractMatrix, B::SwitchingPeriodicArray) = horzcat(SwitchingPeriodicArray(reshape(A,size(A,1),size(A,2),1), [B.dperiod], B.period), B)
hcat(A::SwitchingPeriodicArray, B::AbstractMatrix) = horzcat(A, SwitchingPeriodicArray(reshape(B,size(B,1),size(B,2),1), [A.dperiod], A.period))
hcat(A::AbstractMatrix, B::SwitchingPeriodicArray) = horzcat(SwitchingPeriodicArray(reshape(A,size(A,1),size(A,2),1), [B.dperiod], B.period), B)

# +(A::SwitchingPeriodicArray, C::AbstractMatrix) = +(A, SwitchingPeriodicArray(reshape(C,size(A,1),size(A,2),1), [A.dperiod], A.period))
# +(A::AbstractMatrix, C::SwitchingPeriodicArray) = +(SwitchingPeriodicArray([A], C.ns, [C.dperiod], C.period), C)
# SwitchingPeriodicArray(M::VecOrMat{T}, period::Real; nperiod::Int = 1) where T = SwitchingPeriodicArray(reshape(M,size(M,1),size(M,2),1), [1], period; nperiod)



function vertcat(A::SwitchingPeriodicArray, B::SwitchingPeriodicArray)
    isconstant(A) || isconstant(B) || A.Ts ≈ B.Ts || error("A and B must have the same sampling time")
    period = promote_period(A, B)
    T = promote_type(eltype(A),eltype(B))
    ma, na, pa = size(A.M)
    mb, nb, pb = size(B.M)
    na == nb || throw(DimensionMismatch("number of columnss of A $na not equal to the number of columns of B $nb"))

    if A.period == B.period && A.nperiod == B.nperiod && A.ns == B.ns 
        N = length(A.ns)
        X = Array{T,3}(undef, ma+mb, na, N)
        i1 = 1:ma; i2 = ma+1:ma+mb
        for i = 1:N
            copyto!(view(X,i1,:,i),view(A.M,:,:,i))
            copyto!(view(X,i2,:,i),view(B.M,:,:,i))
        end
        return SwitchingPeriodicArray{:d,T}(X, A.ns, A.period, A.nperiod)
    end
    if isconstant(A) 
        N = length(B.ns)
        X = Array{T,3}(undef, ma+mb, na, N)
        i1 = 1:ma; i2 = ma+1:ma+mb
        for i = 1:N
            copyto!(view(X,i1,:,i),view(A.M,:,:,1))
            copyto!(view(X,i2,:,i),view(B.M,:,:,i))
        end
        return SwitchingPeriodicArray{:d,T}(X, B.ns, B.period, B.nperiod)
    elseif isconstant(B) 
        N = length(A.ns)
        X = Array{T,3}(undef, ma+mb, na, N)
        i1 = 1:ma; i2 = ma+1:ma+mb
        for i = 1:N
            copyto!(view(X,i1,:,i),view(A.M,:,:,i))
            copyto!(view(X,i2,:,i),view(B.M,:,:,1))
        end
        return SwitchingPeriodicArray{:d,T}(X, A.ns, A.period, A.nperiod)
    elseif A.period == B.period 
        nperiod = A.nperiod
        if nperiod == B.nperiod
           ns = unique(sort([A.ns;B.ns]))
        else
           nperiod = gcd(A.nperiod,B.nperiod)
           ns = unique(sort([vcat([(i-1)*A.dperiod .+ A.ns for i in 1:div(A.nperiod,nperiod)]...);
                             vcat([(i-1)*B.dperiod .+ B.ns for i in 1:div(B.nperiod,nperiod)]...)]))
        end
        N = length(ns)   
        X = Array{T,3}(undef, ma+mb, na, N)
        i1 = 1:ma; i2 = ma+1:ma+mb
        for i = 1:N
            copyto!(view(X,i1,:,i),kpmeval(A,ns[i]))
            copyto!(view(X,i2,:,i),kpmeval(B,ns[i]))
        end
        return SwitchingPeriodicArray{:d,promote_type(eltype(A),eltype(B))}(X, ns, A.period, nperiod)   
    else
        # TODO: implement explicit computations
        period = promote_period(A,B)
        return vertcat(set_period(A,period),set_period(B,period))
        #error("periods must be equal for vertical concatenation")
    end
end
Base.vcat(A::SwitchingPeriodicArray, B::SwitchingPeriodicArray) = vertcat(A,B)
vertcat(A::SwitchingPeriodicArray, B::AbstractMatrix) = vertcat(A, SwitchingPeriodicArray(reshape(B,size(B,1),size(B,2),1), [A.dperiod], A.period))
vertcat(A::AbstractMatrix, B::SwitchingPeriodicArray) = vertcat(SwitchingPeriodicArray(reshape(A,size(A,1),size(A,2),1), [B.dperiod], B.period), B)
vcat(A::SwitchingPeriodicArray, B::AbstractMatrix) = vertcat(A, SwitchingPeriodicArray(reshape(B,size(B,1),size(B,2),1), [A.dperiod], A.period))
vcat(A::AbstractMatrix, B::SwitchingPeriodicArray) = vertcat(SwitchingPeriodicArray(reshape(A,size(A,1),size(A,2),1), [B.dperiod], B.period), B)

function blockdiag(A::SwitchingPeriodicArray, B::SwitchingPeriodicArray)
    isconstant(A) || isconstant(B) || A.Ts ≈ B.Ts || error("A and B must have the same sampling time")
    period = promote_period(A, B)
    T = promote_type(eltype(A),eltype(B))
    ma, na, pa = size(A.M)
    mb, nb, pb = size(B.M)

    if A.period == B.period && A.nperiod == B.nperiod && A.ns == B.ns 
        N = length(A.ns)
        X = zeros(T, ma+mb, na+nb, N)
        i1 = 1:ma; i2 = ma+1:ma+mb
        j1 = 1:na; j2 = na+1:na+nb
        for i = 1:N
            copyto!(view(X,i1,j1,i),view(A.M,:,:,i))
            copyto!(view(X,i2,j2,i),view(B.M,:,:,i))
        end
        return SwitchingPeriodicArray{:d,T}(X, A.ns, A.period, A.nperiod)
    end
    if isconstant(A) 
        N = length(B.ns)
        X = zeros(T, ma+mb, na+nb, N)
        i1 = 1:ma; i2 = ma+1:ma+mb
        j1 = 1:na; j2 = na+1:na+nb
        for i = 1:N
            copyto!(view(X,i1,j1,i),view(A.M,:,:,1))
            copyto!(view(X,i2,j2,i),view(B.M,:,:,i))
        end
        return SwitchingPeriodicArray{:d,T}(X, B.ns, B.period, B.nperiod)
    elseif isconstant(B) 
        N = length(A.ns)
        X = zeros(T, ma+mb, na+nb, N)
        i1 = 1:ma; i2 = ma+1:ma+mb
        j1 = 1:na; j2 = na+1:na+nb
        for i = 1:N
            copyto!(view(X,i1,j1,i),view(A.M,:,:,i))
            copyto!(view(X,i2,j2,i),view(B.M,:,:,1))
        end
        return SwitchingPeriodicArray{:d,T}(X, A.ns, A.period, A.nperiod)
    elseif A.period == B.period 
        nperiod = A.nperiod
        if nperiod == B.nperiod
           ns = unique(sort([A.ns;B.ns]))
        else
           nperiod = gcd(A.nperiod,B.nperiod)
           ns = unique(sort([vcat([(i-1)*A.dperiod .+ A.ns for i in 1:div(A.nperiod,nperiod)]...);
                             vcat([(i-1)*B.dperiod .+ B.ns for i in 1:div(B.nperiod,nperiod)]...)]))
        end
        N = length(ns)   
        X = zeros(T, ma+mb, na+nb, N)
        i1 = 1:ma; i2 = ma+1:ma+mb
        j1 = 1:na; j2 = na+1:na+nb
        for i = 1:N
            copyto!(view(X,i1,j1,i),kpmeval(A,ns[i]))
            copyto!(view(X,i2,j2,i),kpmeval(B,ns[i]))
        end
        return SwitchingPeriodicArray{:d,promote_type(eltype(A),eltype(B))}(X, ns, A.period, nperiod)   
    else
        # TODO: implement explicit computations
        period = promote_period(A,B)
        return blockdiag(set_period(A,period),set_period(B,period))
        #error("periods must be equal for bloc diagonal appending")
    end
end


# Operations with periodic function matrices
"""
    pmderiv(A::PeriodicFunctionMatrix; h:Union{Real,Missing} = missing, method = "cd", discont = false) 
    pmderiv(A::PeriodicTimeSeriesMatrix; h:Union{Real,Missing} = missing, method = "cd", discont = false) 
    pmderiv(A::PeriodicSwitchingMatrix; h:Union{Real,Missing} = missing, method = "cd", discont = false) 

Compute the derivative of a continuous-time periodic matrix using finite difference formulas. By default `method = "cd"`, in which case 
the central difference formula is used for approximating the derivative. 
If `method = "4d"`, a fourth order finite difference formula is used for higher accuracy. If `method = ""`, the forward difference formula is used if the 
time difference `h` is missing or `h > 0`, and a backward difference formula is used if `h < 0`. If `discont = true`, then initial discountinuities at `t = 0` and 
terminal discountinuities at `tsub := A.period/A.nperiod` (the subperiod) are 
avoided by using the forward or backward differentiation formulas at `t =  0` and at `t = tsub`, respectively. 
This approach generally leads to lower accuracy estimations at `t = 0` and `t = tsub`.   

_Note:_ To allow the application of the finite difference approximations to periodic matrices of types `PeriodicTimeSeriesMatrix` and `PeriodicSwitchingMatrix`, 
these are converted to `PeriodicFunctionMatrix` type. Due to inherent discountinuities of these representations, the accuracy of derivate estimations is usualy poor.
To increase the accuracy, it is recommended to perform these conversions before calling the `pmderiv` function, by employing splines based interpolation formulas, 
as provided by the function [`ts2pfm`](@ref).  
"""
function pmderiv(A::PeriodicFunctionMatrix{:c,T};  h::Union{Missing,Real} = missing,  method = "cd", discont = false) where {T}
    isconstant(A) && (return PeriodicFunctionMatrix{:c,T}(t -> zeros(T,A.dims...), A.period, A.dims, A.nperiod, true))
    # centered differences
    tsub = A.period/A.nperiod
    # try to avoid discontinuities at the end of the interval
    if method == "cd" 
       ismissing(h) ? (h = tsub*sqrt(eps(float(T)))) : h = abs(h)
       if discont 
          return PeriodicFunctionMatrix{:c,T}(t -> t+h >= tsub ? (tpmeval(A,t)-tpmeval(A,t-h/4))/(h/4) : 
              (t-h < 0 ? (tpmeval(A,t+h/4)-tpmeval(A,t))/(h/4) : (tpmeval(A,t+h)-tpmeval(A,t-h))/(2*h)), A.period, A.dims, A.nperiod, false)
       else
          return PeriodicFunctionMatrix{:c,T}(t -> (tpmeval(A,t+h)-tpmeval(A,t-h))/(2*h), A.period, A.dims, A.nperiod, false)
       end
    end
    # fourth-order differences
    if method == "4d" 
        ismissing(h) ? (h = tsub*sqrt(sqrt(eps(float(T))))) : h = abs(h)
        if discont 
           return PeriodicFunctionMatrix{:c,T}(t -> t+h >= tsub ? (tpmeval(A,t)-tpmeval(A,t-h/4))/(h/4) : 
                  (t-h < 0 ? (tpmeval(A,t+h/4)-tpmeval(A,t))/(h/4) : (-tpmeval(A,t+2h) + tpmeval(A,t-2h)  + 8*(tpmeval(A,t+h) - tpmeval(A,t-h)))/(12h)), A.period, A.dims, A.nperiod, false)
        else
           return PeriodicFunctionMatrix{:c,T}(t -> (-tpmeval(A,t+2h) + tpmeval(A,t-2h)  + 8*(tpmeval(A,t+h) - tpmeval(A,t-h)))/(12h), A.period, A.dims, A.nperiod, false)
        end
    end
    # first-order differences 
    ismissing(h) && (h = tsub*sqrt(eps(float(T))))
    if discont 
       if h > 0
          return PeriodicFunctionMatrix{:c,T}(t -> t+h >= tsub ? (tpmeval(A,t)-tpmeval(A,t-h))/h : (tpmeval(A,t+h)-tpmeval(A,t))/h, A.period, A.dims, A.nperiod, false)
       else
          return PeriodicFunctionMatrix{:c,T}(t -> t+h < 0 ? (tpmeval(A,t)-tpmeval(A,t-h))/h : (tpmeval(A,t+h)-tpmeval(A,t))/h, A.period, A.dims, A.nperiod, false)
       end
    else
       return PeriodicFunctionMatrix{:c,T}(t -> (tpmeval(A,t+h)-tpmeval(A,t))/h, A.period, A.dims, A.nperiod, false)
    end
end
function LinearAlgebra.inv(A::PeriodicFunctionMatrix{:c,T})  where {T}
    return PeriodicFunctionMatrix{:c,T}(t -> inv(tpmeval(A,t)), A.period, A.dims, A.nperiod, A._isconstant)
end
function LinearAlgebra.tr(A::PeriodicFunctionMatrix{:c,T})  where {T}
    return PeriodicFunctionMatrix{:c,T}(t -> [tr(tpmeval(A,t))], A.period, (1,1), A.nperiod, A._isconstant)
end
"""
    trace(A; rtol=sqrt(eps)) -> Atrace

Compute the trace of a continuous-time periodic matrix over one period. 
For a continuous-time periodic matrix `A(t)`, the resulting `Atrace` is the mean value of the integral of the point-wise trace over a complete period. 
The involved time integral are evaluated using the adaptive Gauss-Kronrod quadrature with a relative error tolerance `rtol`. 
"""
function trace(A::Union{HarmonicArray,PeriodicFunctionMatrix}; rtol=sqrt(eps())) 
    isconstant(A) && (return tr(tpmeval(A, 0)))
    tsub = A.period/A.nperiod
    tt, = quadgk(t -> tr(tpmeval(A,t)), 0., tsub; rtol)
    return tt/tsub
end

function LinearAlgebra.transpose(A::PeriodicFunctionMatrix{:c,T})  where {T}
    return PeriodicFunctionMatrix{:c,T}(t -> transpose(tpmeval(A,t)), A.period, (A.dims[2],A.dims[1]), A.nperiod, A._isconstant)
end
function LinearAlgebra.adjoint(A::PeriodicFunctionMatrix{:c,T})  where {T}
    return PeriodicFunctionMatrix{:c,T}(t -> adjoint(tpmeval(A,t)), A.period, (A.dims[2],A.dims[1]), A.nperiod, A._isconstant)
end
function LinearAlgebra.opnorm(A::PeriodicFunctionMatrix, p::Union{Real, Missing} = missing) 
    if ismissing(p)
       return PeriodicFunctionMatrix{:c,eltype(A)}(t -> [norm(tpmeval(A,t))], A.period, (1,1), A.nperiod, A._isconstant)
    else
       return PeriodicFunctionMatrix{:c,eltype(A)}(t -> [opnorm(tpmeval(A,t),p)], A.period, (1,1), A.nperiod, A._isconstant)
    end
end
"""
    norm(A, p::Real=2) 
    norm(A::PeriodicFunctionMatrix, p::Real=2; rtol=sqrt(eps()), atol = 1000*eps()) 

Compute the `p`-norm of the time-dependent Frobenius-norm of a continuous-time periodic matrix over one period. 
For a continuous-time periodic matrix `A(t)`, the resulting `Anorm` is the `p`-norm of the time-varying Frobenius norm of `A(t)`. 
The involved time integrals are evaluated using the adaptive Gauss-Kronrod quadrature with a relative error tolerance `rtol`
and an absolute tolerance `atol`. For `p = Inf`, the computation involves the minimization of the Frobenius norm of `A(t)` using Brent's method.   

_Note:_ For periodic matrices of `PeriodicTimeSeriesMatrix` and `PeriodicSwitchingMatrix` types, 
the existing implicit time griding is employed to evaluate the involved time integrals using the rectangle method. 
"""
function LinearAlgebra.norm(A::Union{HarmonicArray,PeriodicFunctionMatrix}, p::Real = 2; rtol = sqrt(eps()), atol = 1000. *eps()) 
    isconstant(A) && (return norm(tpmeval(A,0),p))
    tsub = A.period/A.nperiod
    if p == 2
       nrm, = quadgk(t -> norm(tpmeval(A,t))^2, 0., tsub; rtol, atol)
       return sqrt(nrm*A.nperiod)
    elseif isinf(p)
        return -optimize(t->-norm(tpmeval(A,t)),0.,tsub,Optim.Brent(),rel_tol = rtol, abs_tol = atol).minimum
    elseif p == 1    
        nrm, = quadgk(t -> norm(tpmeval(A,t)), 0., tsub; rtol, atol)
        return nrm*A.nperiod
    else
        throw(ArgumentError("only p-norms for p = 1, 2, or Inf are supported"))
    end
end
function +(A::PeriodicFunctionMatrix, B::PeriodicFunctionMatrix)
    A.dims == B.dims || throw(DimensionMismatch("A and B must have the same dimensions"))
    # nta = numerator(rationalize(period/A.period))
    # ntb = numerator(rationalize(period/B.period))
    # nperiod = gcd(nta*A.nperiod,ntb*B.nperiod)
    T = promote_type(eltype(A),eltype(B))
    period, nperiod = promote_period2(A, B)
    if isconstant(A) && isconstant(B)
       return PeriodicFunctionMatrix{:c,T}(t -> A.f(0)+B.f(0), period, A.dims, nperiod, true)
    else
       return PeriodicFunctionMatrix{:c,T}(t -> tpmeval(A,t)+tpmeval(B,t), period, A.dims, nperiod, false)
    end
end
+(A::PeriodicFunctionMatrix, C::AbstractMatrix) = +(A, PeriodicFunctionMatrix(C, A.period))
+(A::AbstractMatrix, C::PeriodicFunctionMatrix) = +(PeriodicFunctionMatrix(A, C.period), C)
-(A::PeriodicFunctionMatrix) = PeriodicFunctionMatrix{:c,eltype(A)}(t -> -A.f(t), A.period, A.dims, A.nperiod,A._isconstant)
-(A::PeriodicFunctionMatrix, B::PeriodicFunctionMatrix) = +(A,-B)
-(A::PeriodicFunctionMatrix, C::AbstractMatrix) = +(A,-C)
-(A::AbstractMatrix, C::PeriodicFunctionMatrix) = +(A, -C)
function (+)(A::PeriodicFunctionMatrix, J::UniformScaling{<:Real})
    A.dims[1] == A.dims[2] || throw(DimensionMismatch("matrix is not square: dimensions are $(A.dims)"))
    PeriodicFunctionMatrix{:c,eltype(A)}(t -> tpmeval(A,t)+J, A.period, A.dims, A.nperiod,A._isconstant)
end
(+)(J::UniformScaling{<:Real}, A::PeriodicFunctionMatrix) = +(A,J)
(-)(A::PeriodicFunctionMatrix, J::UniformScaling{<:Real}) = +(A,-J)
(-)(J::UniformScaling{<:Real}, A::PeriodicFunctionMatrix) = +(-A,J)

function *(A::PeriodicFunctionMatrix, B::PeriodicFunctionMatrix)
    A.dims[2] == B.dims[1] || throw(DimensionMismatch("A and B have incompatible dimensions"))
    # nta = numerator(rationalize(period/A.period))
    # ntb = numerator(rationalize(period/B.period))
    # nperiod = gcd(nta*A.nperiod,ntb*B.nperiod)
    T = promote_type(eltype(A),eltype(B))
    period, nperiod = promote_period2(A, B)
    if isconstant(A) && isconstant(B)
       return PeriodicFunctionMatrix{:c,T}(t -> A.f(0)*B.f(0), period, (A.dims[1],B.dims[2]), nperiod, true)
    else
       return PeriodicFunctionMatrix{:c,T}(t -> tpmeval(A,t)*tpmeval(B,t), period, (A.dims[1],B.dims[2]), nperiod, false)
    end
end
*(A::PeriodicFunctionMatrix, C::AbstractMatrix) = *(A, PeriodicFunctionMatrix(C, A.period))
*(A::AbstractMatrix, C::PeriodicFunctionMatrix) = *(PeriodicFunctionMatrix(A, C.period), C)
function *(A::PeriodicFunctionMatrix, C::Real)
    if iszero(C)
       T = eltype(A)
       PeriodicFunctionMatrix{:c,T}(t -> zeros(T,A.dims...), A.period, A.dims, A.nperiod,true) 
    else
       PeriodicFunctionMatrix{:c,eltype(A)}(t -> C.*tpmeval(A,t), A.period, A.dims, A.nperiod,A._isconstant)
    end
end
*(A::Real, C::PeriodicFunctionMatrix) = *(C,A)
/(A::PeriodicFunctionMatrix, C::Real) = *(A, 1/C)
*(J::UniformScaling{<:Real}, A::PeriodicFunctionMatrix) = J.λ*A
*(A::PeriodicFunctionMatrix, J::UniformScaling{<:Real}) = A*J.λ
for (PMF, MF) in ((:pmmuladdsym, :muladdsym!), (:pmmultraddsym, :multraddsym!), (:pmmuladdtrsym,:muladdtrsym!) )
    @eval begin
        function  $PMF(A::PeriodicFunctionMatrix,B::PeriodicFunctionMatrix,C::PeriodicFunctionMatrix, (α,β) = (true, true))
            A.dims[1] == A.dims[2] || throw(ArgumentError("A must be a square matrix"))
            if $MF == muladdsym!
                A.dims[1] == B.dims[1] || throw(ArgumentError("A and B must have the same number of rows"))
                A.dims[1] == C.dims[2] || throw(ArgumentError("A and C must have the same number of columns"))
                B.dims[2] == C.dims[1] || throw(ArgumentError("the number of columns of B must be equal to the number of rows of C"))
            elseif $MF == multraddsym!
                A.dims[1] == B.dims[2] || throw(ArgumentError("the number of rows of A must be equal to the number of columns of B"))
                A.dims[1] == C.dims[2] || throw(ArgumentError("A and C must have the same number of columns"))
                B.dims[1] == C.dims[1] || throw(ArgumentError("B and C must have the same number of rows"))
            else
                A.dims[1] == B.dims[1] || throw(ArgumentError("A and B must have the same number of rows"))
                A.dims[1] == C.dims[1] || throw(ArgumentError("the number of columns of A must be equal with the number of rows of C"))
                B.dims[2] == C.dims[2] || throw(ArgumentError("B and C must have the same number of columns"))
            end
            # nperiod = gcd(A.nperiod,B.nperiod,C.nperiod)
            T = promote_type(eltype(A),eltype(B),eltype(C))
            period, nperiod = promote_period2(A, B, C)
            if isconstant(A) && isconstant(B) && isconstant(C)
               return PeriodicFunctionMatrix{:c,T}(t -> $MF(T.(copy(A.f(0))),B.f(0),C.f(0),(α,β)), period, (A.dims[1],A.dims[2]), nperiod, true)
            else
               return PeriodicFunctionMatrix{:c,T}(t -> $MF(T.(copy(tpmeval(A,t))),tpmeval(B,t),tpmeval(C,t),(α,β)), period, (A.dims[1],A.dims[2]), nperiod, false)
            end
        end
        $PMF(A::PeriodicFunctionMatrix,B::PeriodicFunctionMatrix,C::PeriodicFunctionMatrix, α, β) = $PMF(A, B, C, (α,β))
    end
end
# function muladdtrsym!(A::AbstractMatrix, B::AbstractMatrix, C::AbstractMatrix, (α,β) = (true, true))
#     # compute in A the symmetrix matrix α*A +  β*transpose(B)*C
#     n = LinearAlgebra.checksquare(A)
#     n == size(B,1) || throw(ArgumentError("A and B must have the same number of rows"))
#     n == size(C,1) || throw(ArgumentError("A and C must have the same number of rows"))
#     m = size(B,2)
#     m == size(C,2) || throw(ArgumentError("B and C have incompatible dimensions"))
#     ZERO = zero(promote_type(eltype(B),eltype(C)))
#     if α == 0
#         for i = 1:n
#             for j = i:n
#                 temp = ZERO
#                 for k = 1:m
#                     temp += (B[i,k]*C[j,k])
#                 end
#                 A[i,j] = β*temp
#                 A[j,i] = A[i,j]
#             end
#         end
#     else
#         for i = 1:n
#             for j = i:n
#                 temp = ZERO
#                 for k = 1:m
#                     temp += (B[i,k]*C[j,k])
#                 end
#                 A[i,j] = α*A[i,j]+β*temp
#                 A[j,i] = A[i,j]
#             end
#         end
#     end
#     return A
# end
function mulsym(B::AbstractMatrix, C::AbstractMatrix, β = true)
    # compute the symmetrix matrix A = β*B*C
    n, m = size(B) 
    (m, n) == size(C) || throw(ArgumentError("B' and C must have the same dimensions"))
    T = promote_type(eltype(B),eltype(C))
    A = Matrix{T}(undef,n,n)
    ZERO = zero(T)
    for i = 1:n
        for j = i:n
            temp = ZERO
            for k = 1:m
                temp += (B[i,k]*C[k,j])
            end
            A[i,j] = β*temp
            A[j,i] = A[i,j]
        end
    end
    return A
end

function multrsym(B::AbstractMatrix, C::AbstractMatrix, β = true)
    # compute the symmetrix matrix A = β*B*transpose(C)
    n, m = size(B) 
    (n, m) == size(C) || throw(ArgumentError("B and C must have the same dimensions"))
    T = promote_type(eltype(B),eltype(C))
    A = Matrix{T}(undef,n,n)
    ZERO = zero(T)
    for i = 1:n
        for j = i:n
            temp = ZERO
            for k = 1:m
                temp += (B[i,k]*C[j,k])
            end
            A[i,j] = β*temp
            A[j,i] = A[i,j]
        end
    end
    return A
end
function trmulsym(B::AbstractMatrix, C::AbstractMatrix, β = true)
    # compute the symmetrix matrix A = β*transpose(B)*C
    m, n = size(B) 
    (m, n) == size(C) || throw(ArgumentError("B and C must have the same dimensions"))
    T = promote_type(eltype(B),eltype(C))
    A = Matrix{T}(undef,n,n)
    ZERO = zero(T)
    for i = 1:n
        for j = i:n
            temp = ZERO
            for k = 1:m
                temp += (B[k,i]*C[k,j])
            end
            A[i,j] = β*temp
            A[j,i] = A[i,j]
        end
    end
    return A
end


for (PMF, MF) in ((:pmmulsym, :mulsym), (:pmtrmulsym, :trmulsym), (:pmmultrsym,:multrsym) )
    @eval begin
        function $PMF(B::PeriodicFunctionMatrix,C::PeriodicFunctionMatrix, β = true)
            # compute the symmetric results βB*C, βB'*C, and βB*C'. 
            # period = promote_period(B, C)
            # nperiod = gcd(B.nperiod,C.nperiod)
            T = promote_type(eltype(B),eltype(C))
            period, nperiod = promote_period2(B, C)
            n = $PMF == pmtrmulsym ? size(B,2) : size(B,1)
            if isconstant(B) && isconstant(C)
               return PeriodicFunctionMatrix{:c,T}(t -> $MF(B.f(0),C.f(0),β), period, (n,n), nperiod, true)
            else
               return PeriodicFunctionMatrix{:c,T}(t -> $MF(tpmeval(B,t),tpmeval(C,t),β), period, (n,n), nperiod, false)
            end
        end
    end
end

for PMF in (:pmmuladdsym, :pmmultraddsym, :pmmuladdtrsym)
    for PM in (:HarmonicArray,)
        @eval begin
            $PMF(A::$PM,B::$PM,C::$PM, (α,β) = (true, true)) = convert($PM,$PMF(convert(PeriodicFunctionMatrix,A), convert(PeriodicFunctionMatrix,B), convert(PeriodicFunctionMatrix,C), (α,β)))
            $PMF(A::$PM,B::$PM,C::$PM, α, β) = $PMF(A, B, C, (α,β))
            # $PMF(A::$PM,B::AbstractMatrix,C::$PM, (α,β) = (true, true)) = $PMF(A, $PM(B, A.period), C, (α,β))
            # $PMF(A::$PM,B::$PM,C::AbstractMatrix, (α,β) = (true, true)) = $PMF(A, B, $PM(C, A.period), (α,β))
            # $PMF(A::$PM,B::AbstractMatrix,C::AbstractMatrix, (α,β) = (true, true)) = $PMF(A, $PM(B, A.period), $PM(C, A.period), (α,β))
            # $PMF(A::AbstractMatrix,B::$PM,C::$PM, (α,β) = (true, true)) = $PMF($PM(A, B.period), B, C, (α,β))
            # $PMF(A::AbstractMatrix,B::AbstractMatrix,C::$PM, (α,β) = (true, true)) = $PMF($PM(A, C.period), $PM(B, C.period), C, (α,β))
            # $PMF(A::AbstractMatrix,B::$PM,C::AbstractMatrix, (α,β) = (true, true)) = $PMF($PM(A, B.period), B, $PM(C, B.period), (α,β))
        end
    end
end


for PMF in (:pmmuladdsym, :pmmultraddsym, :pmmuladdtrsym)
    for PM in (:PeriodicFunctionMatrix, :HarmonicArray)
        @eval begin
            $PMF(A::$PM,B::AbstractMatrix,C::$PM, (α,β) = (true, true)) = $PMF(A, $PM(B, A.period), C, (α,β))
            $PMF(A::$PM,B::AbstractMatrix,C::$PM, α, β) = $PMF(A, B, C, (α,β))
            $PMF(A::$PM,B::$PM,C::AbstractMatrix, (α,β) = (true, true)) = $PMF(A, B, $PM(C, A.period), (α,β))
            $PMF(A::$PM,B::$PM,C::AbstractMatrix, α, β) = $PMF(A, B, C, (α,β))
            $PMF(A::$PM,B::AbstractMatrix,C::AbstractMatrix, (α,β) = (true, true)) = $PMF(A, $PM(B, A.period), $PM(C, A.period), (α,β))
            $PMF(A::$PM,B::AbstractMatrix,C::AbstractMatrix, α, β) = $PMF(A, B, C, (α,β))
            $PMF(A::AbstractMatrix,B::$PM,C::$PM, (α,β) = (true, true)) = $PMF($PM(A, B.period), B, C, (α,β))
            $PMF(A::AbstractMatrix,B::$PM,C::$PM, α, β) = $PMF(A, B, C, (α,β))
            $PMF(A::AbstractMatrix,B::AbstractMatrix,C::$PM, (α,β) = (true, true)) = $PMF($PM(A, C.period), $PM(B, C.period), C, (α,β))
            $PMF(A::AbstractMatrix,B::AbstractMatrix,C::$PM, α, β) = $PMF(A, B, C, (α,β))
            $PMF(A::AbstractMatrix,B::$PM,C::AbstractMatrix, (α,β) = (true, true)) = $PMF($PM(A, B.period), B, $PM(C, B.period), (α,β))
            $PMF(A::AbstractMatrix,B::$PM,C::AbstractMatrix, α, β) = $PMF(A, B, C, (α,β))
            # $PMF(A::$PM,B::AbstractMatrix,C::$PM, (α,β) = (true, true)) = $PMF(A, $PM(B, A.period), C, (α,β))
            # $PMF(A::$PM,B::$PM,C::AbstractMatrix, (α,β) = (true, true)) = $PMF(A, B, $PM(C, A.period), (α,β))
            # $PMF(A::$PM,B::AbstractMatrix,C::AbstractMatrix, (α,β) = (true, true)) = $PMF(A, $PM(B, A.period), $PM(C, A.period), (α,β))
            # $PMF(A::AbstractMatrix,B::$PM,C::$PM, (α,β) = (true, true)) = $PMF($PM(A, B.period), B, C, (α,β))
            # $PMF(A::AbstractMatrix,B::AbstractMatrix,C::$PM, (α,β) = (true, true)) = $PMF($PM(A, C.period), $PM(B, C.period), C, (α,β))
            # $PMF(A::AbstractMatrix,B::$PM,C::AbstractMatrix, (α,β) = (true, true)) = $PMF($PM(A, B.period), B, $PM(C, B.period), (α,β))
        end
    end
end
for PMF in (:pmmulsym, :pmtrmulsym, :pmmultrsym)
    for PM in (:HarmonicArray,:PeriodicTimeSeriesMatrix)
        @eval begin
            $PMF(B::$PM,C::$PM, β = true) = convert($PM,$PMF(convert(PeriodicFunctionMatrix,B), convert(PeriodicFunctionMatrix,C), β))
        end
    end
end


function horzcat(A::PeriodicFunctionMatrix, B::PeriodicFunctionMatrix)
    A.dims[1] == B.dims[1] || throw(DimensionMismatch("A and B have incompatible row dimensions"))
    # nperiod = numerator(rationalize(period/A.period))*A.nperiod
    T = promote_type(eltype(A),eltype(B))
    period, nperiod = promote_period2(A, B)
    if isconstant(A) && isconstant(B)
       return PeriodicFunctionMatrix{:c,T}(t -> [A.f(0) B.f(0)], period, (A.dims[1],A.dims[2]+B.dims[2]), nperiod, true)
    else
       return PeriodicFunctionMatrix{:c,T}(t -> [tpmeval(A,t) tpmeval(B,t)], period, (A.dims[1],A.dims[2]+B.dims[2]), nperiod, false)
    end
end
horzcat(A::PeriodicFunctionMatrix, C::AbstractMatrix) = horzcat(A, PeriodicFunctionMatrix(C, A.period))
horzcat(A::AbstractMatrix, C::PeriodicFunctionMatrix) = horzcat(PeriodicFunctionMatrix(A, C.period), C)
Base.hcat(A::PeriodicFunctionMatrix, B::PeriodicFunctionMatrix) = horzcat(A,B)
Base.hcat(A::PeriodicFunctionMatrix, B::AbstractMatrix) = horzcat(A,B)
Base.hcat(A::AbstractMatrix, B::PeriodicFunctionMatrix) = horzcat(A,B)

function vertcat(A::PeriodicFunctionMatrix, B::PeriodicFunctionMatrix)
    #period = promote_period(A, B)
    A.dims[2] == B.dims[2] || throw(DimensionMismatch("A and B have incompatible column dimensions"))
    #nperiod = numerator(rationalize(period/A.period))*A.nperiod
    T = promote_type(eltype(A),eltype(B))
    period, nperiod = promote_period2(A, B)
    if isconstant(A) && isconstant(B)
       return PeriodicFunctionMatrix{:c,T}(t -> [A.f(0); B.f(0)], period, (A.dims[1]+B.dims[1],A.dims[2]), nperiod, true)
    else
       return PeriodicFunctionMatrix{:c,T}(t -> [tpmeval(A,t); tpmeval(B,t)], period, (A.dims[1]+B.dims[1],A.dims[2]), nperiod, false)
    end
end
Base.vcat(A::PeriodicFunctionMatrix, B::PeriodicFunctionMatrix) = vertcat(A,B)
vertcat(A::PeriodicFunctionMatrix, C::AbstractMatrix) = vertcat(A, PeriodicFunctionMatrix(C, A.period))
vertcat(A::AbstractMatrix, C::PeriodicFunctionMatrix) = vertcat(PeriodicFunctionMatrix(A, C.period), C)
Base.vcat(A::PeriodicFunctionMatrix, B::AbstractMatrix) = vertcat(A,B)
Base.vcat(A::AbstractMatrix, B::PeriodicFunctionMatrix) = vertcat(A,B)

function blockdiag(A::PeriodicFunctionMatrix, B::PeriodicFunctionMatrix)
    #period = promote_period(A, B)
    #nperiod = numerator(rationalize(period/A.period))*A.nperiod
    T = promote_type(eltype(A),eltype(B))
    period, nperiod = promote_period2(A, B)
    if isconstant(A) && isconstant(B)
       return PeriodicFunctionMatrix{:c,T}(t -> bldiag(A.f(0), B.f(0)), period, (A.dims[1]+B.dims[1],A.dims[2]+B.dims[2]), nperiod, true)
    else
       return PeriodicFunctionMatrix{:c,T}(t -> bldiag(tpmeval(A,t),tpmeval(B,t)), period, (A.dims[1]+B.dims[1],A.dims[2]+B.dims[2]), nperiod, false)
    end
end




#Base.iszero(A::PeriodicFunctionMatrix) = iszero(A.f(rand()*A.period))
function PeriodicMatrices.iszero(f::Function, period::Real; rtol::Float64 = eps(), atol::Float64 = sqrt(eps()))  
    return abs(optimize(t->-norm(f(t),Inf),0,period,Optim.Brent(),rel_tol = rtol).minimum) < atol 
end
 
function Base.iszero(A::PeriodicFunctionMatrix; rtol::Float64 = eps(), atol::Float64 = sqrt(eps()))  
    #return abs(optimize(t->-norm(A.f(t),Inf),0,A.period/A.nperiod,Optim.Brent(),rel_tol = rtol).minimum) < atol 
    return iszero(A.f, A.period/A.nperiod; rtol, atol)
end
 
LinearAlgebra.issymmetric(A::PeriodicFunctionMatrix) = issymmetric(A.f(rand()*A.period))
function ==(A::PeriodicFunctionMatrix, B::PeriodicFunctionMatrix)
    (A.period == B.period && A.nperiod == B.nperiod) || (return false) 
    isconstant(A) && isconstant(B) && (return isequal(A.f(0), B.f(0)))
    ts = rand()*A.period/A.nperiod
    isequal(A.f(ts), B.f(ts)) && iszero(A-B)
end
function Base.isapprox(A::PeriodicFunctionMatrix, B::PeriodicFunctionMatrix; kwargs...)
    (A.period == B.period && A.nperiod == B.nperiod) || (return false) 
    isconstant(A) && isconstant(B) && (return isapprox(A.f(0), B.f(0); kwargs...))
    tsi = rand(10)*A.period/A.nperiod
    for ts in tsi
        isapprox(A.f(ts), B.f(ts); kwargs...)  || (return false)
    end
    return true
end
function Base.isapprox(A::PeriodicFunctionMatrix, J::UniformScaling{<:Real}; atol::Real=0., rtol::Real=atol>0 ? 0. : sqrt(eps(eltype(A))))
    return isconstant(A;check_const=true,atol, rtol) && isapprox(A.f(0), J; atol, rtol)
end
Base.isapprox(J::UniformScaling{<:Real}, A::PeriodicFunctionMatrix; kwargs...) = isapprox(A, J; kwargs...)

# Operations with harmonic arrays
# function pmrand(::Type{T}, n::Int, m::Int, period::Real = 2*pi; nh::Int = 1) where {T}
#     HarmonicArray(rand(T,n,m), [rand(T,n,m) for i in 1:nh], [rand(T,n,m) for i in 1:nh], period) 
# end    
"""
    pmrand(::Type{PM}, n, m[, period = 2pi]; nh = 1) 
    pmrand(::Type{PM{:c,T}}, n, m[, period = 2pi]; nh = 1) 
    pmrand(n, m[, period = 2pi]; nh = 1) 
    pmrand(PeriodicTimeSeriesMatrix, n, m[, period = 2pi]; ns = 1) 
    pmrand(PeriodicTimeSeriesMatrix{:c,T}, n, m[, period = 2pi]; ns = 1) 
    pmrand(PeriodicSwitchingMatrix, n, m[, period = 2pi]; ts = [0]) 
    pmrand(PeriodicSwitchingMatrix{:c,T}, n, m[, period = 2pi]; ts = [0]) 

Generate a random `n×m` continuous-time periodic matrix of type `PM ` or `PM{:c,T} ` with period  `period ` (default:  `period = 2pi `).
`PM` specifies the resulting type of the generated periodic matrix. For `PM = HarmonicArray`, or `PM = PeriodicFunctionMatrix`, or `PM = PeriodicSymbolicMatrix`
or `PM = FourierFunctionMatrix`, the resulting periodic matrix corresponds 
to a random harmonic representation with  `nh ` harmonic components (default:  `nh = 1 `). 
The type  `T` of matrix elements can be specified using, e.g. `HarmonicArray{:c,T}` instead `HarmonicArray`, 
which assumes by default `T = Float64`. If `PM` is omitted, then by default `PM = HarmonicArray`.

For `PM = PeriodicTimeSeriesMatrix`, `ns` specifies the number of component matrices.

For `PM = PeriodicSwitchingMatrix`, the vector `ts` specifies the switching times. 
"""
function pmrand(n::Int, m::Int, period::Real = 2*pi; nh::Int = 1)
    pmrand(HarmonicArray{:c,Float64}, n, m, period; nh)
end
function pmrand(::Type{PM}, n::Int, m::Int, period::Real = 2*pi; nh::Int = 1) where 
          {T,PM <: Union{HarmonicArray{:c,T},PeriodicFunctionMatrix{:c,T}}}
    A = HarmonicArray(rand(T,n,m), [rand(T,n,m) for i in 1:nh], [rand(T,n,m) for i in 1:nh], period) 
    PM <: HarmonicArray && (return A)
    return convert(PM,A)
end 
function pmrand(::Type{PM}, n::Int, m::Int, period::Real = 2*pi; kwargs...) where 
         {PM <: Union{HarmonicArray,PeriodicFunctionMatrix}}
    pmrand(PM{:c,Float64}, n, m, period; kwargs...)
end 
function pmrand(::Type{PM}, n::Int, m::Int, period::Real = 2*pi; ns::Int = 1) where {T,PM <: PeriodicTimeSeriesMatrix{:c,T}}
    return PeriodicTimeSeriesMatrix{:c,T}([rand(T,n,m) for i in 1:ns], period) 
end 
function pmrand(::Type{PM}, n::Int, m::Int, period::Real = 2*pi; ns::Int = 1) where {PM <: PeriodicTimeSeriesMatrix}
    return PeriodicTimeSeriesMatrix{:c,Float64}([rand(n,m) for i in 1:ns], period) 
end 
function pmrand(::Type{PM}, n::Int, m::Int, period::Real = 2*pi; ts::Vector{<:Real} = [0.]) where {T,PM <: PeriodicSwitchingMatrix{:c,T}}
    return PeriodicSwitchingMatrix{:c,T}([rand(T,n,m) for i in 1:length(ts)], ts, period) 
end 
function pmrand(::Type{PM}, n::Int, m::Int, period::Real = 2*pi; ts::Vector{<:Real} = [0.]) where {PM <: PeriodicSwitchingMatrix}
    return PeriodicSwitchingMatrix{:c,Float64}([rand(n,m) for i in 1:length(ts)], ts, period) 
end 

"""
    pmrand(::Type{PM}, n::Int, m::Int[, period = 10]; ns = 1) 
    pmrand(::Type{PM{:d,T}}, n::Int, m::Int[, period = 10]; ns = 1) 
    pmrand(::Type{PM}, n::Vector{Int}, m::Vector{Int}[, period = 1]) 
    pmrand(SwitchingPeriodicMatrix, n, m[, period = 10]; ns = [period]) 
    pmrand(SwitchingPeriodicArray, n, m[, period = 10]; ns = [period]) 

Generate a random `n×m` discrete-time periodic matrix of type `PM` or `PM{:d,T}` with period  `period ` (default:  `period = 10`)
with  `ns` component matrices (default: `ns = 1`). 
If `PM = PeriodicMatrix` or `PM = PeriodicArray`, `ns` specifies the number of component matrices (default: `ns = 10`).
If `PM = PeriodicMatrix`, then two integer vectors `n` and `m` containing the row and column dimensions of the
the component matrices, respectively, can be used to specify periodic matrices with time-varying dimensions. 

If `PM = SwitchingPeriodicMatrix` or `PM = SwitchingPeriodicArray`, the integer vector `ns` specifies the switching moments (default: `ns = [period]).
The type  `T` of matrix elements can be specified using, e.g. `PeriodicMatrix{:d,T}` instead `PeriodicMatrix`, 
which assumes by default `T = Float64`.
"""
function pmrand(::Type{PM},m::Vector{Int},n::Vector{Int}, period::Real = 10) where {PM <: PeriodicMatrix}
    lm = length(m)
    ln = length(n)
    return PeriodicMatrix{:d,Float64}([rand(m[mod(i-1,lm)+1], n[mod(i-1,ln)+1]) for i in 1:lcm(lm,ln)],period)
end
function pmrand(::Type{PM}, n::Int, m::Int, period::Real = 10; ns::Int = 1) where {T,PM <: Union{PeriodicMatrix{:d,T},PeriodicArray{:d,T}}}
    if PM <: PeriodicMatrix 
        PeriodicMatrix{:d,T}([rand(T,n,m) for i in 1:ns], period) 
    else
        PeriodicArray{:d,T}(rand(T,n,m,ns), period) 
    end
end 
pmrand(::Type{PM}, n::Int, m::Int, period::Real = 10; kwargs...) where {PM <: Union{PeriodicMatrix,PeriodicArray,SwitchingPeriodicMatrix,SwitchingPeriodicArray}} =
    pmrand(PM{:d,Float64}, n, m, period; kwargs...)

function pmrand(::Type{PM}, n::Int, m::Int, period::Real = 10; ns::Vector{Int} = [period]) where {T,PM <: Union{SwitchingPeriodicMatrix{:d,T},SwitchingPeriodicArray{:d,T}}}
    if PM <: SwitchingPeriodicMatrix
        SwitchingPeriodicMatrix{:d,T}([rand(T,n,m) for i in 1:length(ns)], ns, period) 
    else
        SwitchingPeriodicArray{:d,T}(rand(T,n,m,length(ns)), ns, period) 
    end
end 

"""
    pmrand([::Type{T},] m::Vector{Int},n::Vector{Int}[, period = 10])
 
Generate a random discrete-time periodic matrix of type  `PeriodicMatrix` with period  `period ` (default:  `period = 10`). 
The time-varying row and column dimensions of component matrices are specified by the integer vectors `m` and `n`, respectively.
`T` is the type of matrix elements, which assumes by default value `T = Float64` if omitted.
"""
function pmrand(::Type{T},m::Vector{Int},n::Vector{Int}, period::Real = 10) where {T}
    lm = length(m)
    ln = length(n)
    return PeriodicMatrix{:d,T}([rand(T,m[mod(i-1,lm)+1], n[mod(i-1,ln)+1]) for i in 1:lcm(lm,ln)],period)
end
pmrand(m::Vector{Int},n::Vector{Int}, period::Real = 10) = pmrand(Float64,m,n,period)
"""
    pmderiv(A::HarmonicArray) 

Compute the derivative of a continuous-time periodic matrix in harmonic represention.     
"""
function pmderiv(A::HarmonicArray{:c,T}) where {T}
    m, n, l = size(A.values)
    isconstant(A) && (return HarmonicArray{:c,T}(zeros(Complex{T}, m, n, 1), A.period, A.nperiod)) 
    Ahr = similar(A.values)
    ω = 2*pi*A.nperiod/A.period
    Ahr[:,:,1] = zeros(T, m, n)
    for i = 1:l-1
        Ahr[:,:,i+1] .= complex.(imag(A.values[:,:,i+1])*(i*ω),real(A.values[:,:,i+1])*(-i*ω)) 
    end
    return HarmonicArray{:c,T}(Ahr, A.period, nperiod = A.nperiod)
end
function LinearAlgebra.inv(A::HarmonicArray)
    convert(HarmonicArray,inv(convert(PeriodicFunctionMatrix,A)))
    #convert(HarmonicArray,inv(convert(PeriodicTimeSeriesMatrix,A)))
    # ts = (0:N-1)*A.period/A.nperiod/N
    # convert(HarmonicArray,PeriodicTimeSeriesMatrix(inv.(tvmeval(A,collect(ts))),A.period;nperiod = A.nperiod))
    #convert(HarmonicArray,inv(convert(PeriodicFunctionMatrix,A)))
end
function tr(A::HarmonicArray{:c,T}) where {T}
    l = size(A.values,3)
    Ahr = zeros(Complex{T}, 1, 1, l)
    for i = 1:l
        #Ahr[:,:,i] .= complex(tr(real(A.values[:,:,i])),tr(imag(A.values[:,:,i]))) 
        Ahr[1,1,i] = tr(A.values[:,:,i]) 
    end
    return HarmonicArray{:c,T}(Ahr, A.period, nperiod = A.nperiod)
end
# function trace(A::HarmonicArray; K = 128) 
#     isconstant(A) && (return tr(tpmeval(A, 0)))
#     ts = zero(eltype(A))
#     Δ = A.period/A.nperiod/K
#     tt = zero(eltype(Δ))
#     for i = 1:K
#         tt += tr(tpmeval(A, ts))*Δ
#         ts += Δ
#     end 
#     return tt*A.nperiod/A.period
# end

function LinearAlgebra.transpose(A::HarmonicArray{:c,T}) where {T}  
    m, n, l = size(A.values)
    Ahr = similar(Array{Complex{T},3},n,m,l)
    for i = 1:l
        Ahr[:,:,i] .= copy(transpose(A.values[:,:,i])) 
    end
    return HarmonicArray{:c,T}(Ahr, A.period, nperiod = A.nperiod)
end
LinearAlgebra.adjoint(A::HarmonicArray) = transpose(A)
function LinearAlgebra.opnorm(A::HarmonicArray, p::Union{Real, Missing} = missing) 
    return convert(HarmonicArray,opnorm(convert(PeriodicFunctionMatrix,A),p))
end
# function LinearAlgebra.norm(A::HarmonicArray, p::Real = 2; K = 128) 
#     isconstant(A) && (return norm(tpmeval(A, 0)))
#     nrm = zero(eltype(A))
#     Δ = A.period/A.nperiod/K
#     ts = zero(eltype(Δ))
#     if p == 2
#        for i = 1:K
#            nrm += norm(tpmeval(A, ts))^2*Δ
#            ts += Δ
#        end 
#        return sqrt(nrm*A.nperiod)
#     elseif isinf(p)
#         for i = 1:K
#             nrm = max(nrm,norm(tpmeval(A, ts)))
#             ts += Δ
#         end 
#         return nrm
#     elseif p == 1    
#         for i = 1:K
#             nrm += norm(tpmeval(A, ts))*Δ
#             ts += Δ
#         end 
#         return nrm*A.nperiod
#     else
#         throw(ArgumentError("only p-norms for p = 1, 2, or Inf are supported"))
#     end
# end
function +(A::HarmonicArray, B::HarmonicArray)
    if A.period == B.period && A.nperiod == B.nperiod
       m, n, la = size(A.values)
       mb, nb, lb = size(B.values)
       (m, n) == (mb, nb) || throw(DimensionMismatch("A and B must have the same size"))
       T = promote_type(eltype(A.values),eltype(B.values))
       lmax = max(la,lb)
       Ahr = zeros(T,m,n,lmax)
       if la >= lb
          copyto!(view(Ahr,1:m,1:n,1:la),A.values) 
          #Ahr[:,:,1:la] = copy(A.values) 
          Ahr[:,:,1:lb] .+= view(B.values,1:m,1:n,1:lb)  
       else
          copyto!(view(Ahr,1:m,1:n,1:lb),B.values) 
          #Ahr = copy(B.values) 
          Ahr[:,:,1:la] .+= view(A.values,1:m,1:n,1:la)  
       end
       tol = 10*eps(real(T))*max(norm(A.values,Inf),norm(B.values,Inf)) 
       l = lmax
       for i = lmax:-1:2
           norm(Ahr[:,:,i],Inf) > tol && break
           l -= 1
       end
       l < lmax && (Ahr = Ahr[:,:,1:l])
       return HarmonicArray{:c,real(T)}(Ahr, A.period, nperiod = A.nperiod) 
    else
       #TODO: fix different numbers of subperiods 
       convert(HarmonicArray,convert(PeriodicFunctionMatrix,A) + convert(PeriodicFunctionMatrix,B))
    end
end
+(A::HarmonicArray, C::AbstractMatrix) = +(A, HarmonicArray(C, A.period))
+(A::AbstractMatrix, C::HarmonicArray) = +(HarmonicArray(A, C.period), C)
+(A::HarmonicArray, C::PeriodicFunctionMatrix) = +(convert(PeriodicFunctionMatrix,A), C)
+(A::PeriodicFunctionMatrix, C::HarmonicArray) = +(A, convert(PeriodicFunctionMatrix,C))
-(A::HarmonicArray) = HarmonicArray(-A.values, A.period; nperiod = A.nperiod)
-(A::HarmonicArray, B::HarmonicArray) = +(A,-B)
-(A::HarmonicArray, C::AbstractMatrix) = +(A,-C)
-(A::AbstractMatrix, C::HarmonicArray) = +(A, -C)
-(A::HarmonicArray, C::PeriodicFunctionMatrix) = +(A,-C)
-(A::PeriodicFunctionMatrix, C::HarmonicArray) = +(A, -C)
function (+)(A::HarmonicArray, J::UniformScaling{<:Real}) 
    m, n = size(A)
    n == m || throw(DimensionMismatch("matrix is not square: dimensions are $((m,n))"))
    A+Matrix(J(n))
end
(+)(J::UniformScaling{<:Real}, A::HarmonicArray) = +(A,J)
(-)(A::HarmonicArray, J::UniformScaling{<:Real}) = +(A,-J)
(-)(J::UniformScaling{<:Real}, A::HarmonicArray) = +(-A,J)

*(A::HarmonicArray, B::HarmonicArray) = convert(HarmonicArray,convert(PeriodicFunctionMatrix,A) * convert(PeriodicFunctionMatrix,B))


# TODO: perform * explicitly
function *(A::HarmonicArray, C::AbstractMatrix)
    m, n, k = size(A.values)
    nc = size(C,2)
    T = promote_type(eltype(A.values),eltype(C))
    vals = Array{T,3}(undef,m,nc,k)
    [vals[:,:,i] = view(A.values,1:m,1:n,i)*C for i in 1:k]
    return HarmonicArray(vals, A.period; nperiod = A.nperiod)
    #*(A, HarmonicArray(C, A.period))
end
function *(A::AbstractMatrix, C::HarmonicArray) 
    m, n, k = size(C.values)
    ma = size(A,1)
    return HarmonicArray(reshape(A*reshape(C.values,m,n*k),ma,n,k), C.period; nperiod = C.nperiod)
end
*(A::HarmonicArray, C::Real) = HarmonicArray(C*A.values, A.period; nperiod = A.nperiod)
*(C::Real, A::HarmonicArray) = HarmonicArray(C*A.values, A.period; nperiod = A.nperiod)
/(A::HarmonicArray, C::Real) = *(A, 1/C)
*(J::UniformScaling{<:Real}, A::HarmonicArray) = J.λ*A
*(A::HarmonicArray, J::UniformScaling{<:Real}) = A*J.λ
*(A::HarmonicArray, C::PeriodicFunctionMatrix) = *(convert(PeriodicFunctionMatrix,A), C)
*(A::PeriodicFunctionMatrix, C::HarmonicArray) = *(A, convert(PeriodicFunctionMatrix,C))
# function pmmultraddsym(A::HarmonicArray,B::HarmonicArray,C::HarmonicArray, (α,β) = (true, true))
#     convert(HarmonicArray,pmmultraddsym(convert(PeriodicFunctionMatrix,A), convert(PeriodicFunctionMatrix,B), convert(PeriodicFunctionMatrix,C), (α,β)))
# end
# function pmmultraddsym(A::AbstractMatrix,B::HarmonicArray,C::HarmonicArray, (α,β) = (true, true))
#     convert(HarmonicArray,pmmultraddsym(PeriodicFunctionMatrix(A, B.period), convert(PeriodicFunctionMatrix,B), convert(PeriodicFunctionMatrix,C), (α,β)))
# end




Base.iszero(A::HarmonicArray) = all([iszero(A.values[:,:,i]) for i in 1:size(A.values,3)])
LinearAlgebra.issymmetric(A::HarmonicArray) = all([issymmetric(A.values[:,:,i]) for i in 1:size(A.values,3)])
function ==(A::HarmonicArray, B::HarmonicArray)
    (A.period == B.period && A.nperiod == B.nperiod) || (return false) 
    return isequal(A.values, B.values)
end
function Base.isapprox(A::HarmonicArray, B::HarmonicArray; rtol::Real = sqrt(eps(Float64)), atol::Real = 0)
    (A.period == B.period && A.nperiod == B.nperiod) || (return false) 
    isconstant(A) && isconstant(B) && (return isapprox(A.values, B.values; rtol, atol) )
    na = size(A.values,3)
    nb = size(B.values,3)
    if na == nb
       return isapprox(A.values, B.values; rtol, atol) 
    elseif na > nb
        tol = atol+rtol*max(norm(A,1),norm(B,1))
        return all([isapprox(A.values[:,:,i], B.values[:,:,i]; rtol, atol) for i in 1:nb]) && 
               all([norm(A.values[:,:,i],1) < tol for i in nb+1:na])  
    else
        tol = atol+rtol*max(norm(A,1),norm(B,1))
        return all([isapprox(A.values[:,:,i], B.values[:,:,i]; rtol, atol) for i in 1:na]) && 
               all([norm(B.values[:,:,i],1) < tol for i in na+1:nb])  
    end
end
function Base.isapprox(A::HarmonicArray, J::UniformScaling{<:Real}; atol::Real=0., rtol::Real=atol>0 ? 0. : sqrt(eps(eltype(A))))
    return (isconstant(A) || norm(view(A.values,:,:,2:size(A.values,3)),Inf) < atol+rtol) && isapprox(tpmeval(A,0), J; atol, rtol)
end
Base.isapprox(J::UniformScaling{<:Real}, A::HarmonicArray; kwargs...) = isapprox(A, J; kwargs...)

function horzcat(A::HarmonicArray, B::HarmonicArray)
    if A.period == B.period && A.nperiod == B.nperiod
       m, n, la = size(A.values)
       mb, nb, lb = size(B.values)
       m == mb || throw(DimensionMismatch("A and B must have the same number of rows"))
       T = promote_type(eltype(A),eltype(B))
       lmax = max(la,lb)
       Ahr = zeros(Complex{T},m,n+nb,lmax)
       copyto!(view(Ahr,1:m,1:n,1:la), view(A.values,1:m,1:n,1:la))
       copyto!(view(Ahr,1:m,n+1:n+nb,1:lb), view(B.values,1:m,1:nb,1:lb))
       return HarmonicArray{:c,T}(Ahr, A.period, nperiod = A.nperiod) 
    else
       #TODO: fix different numbers of subperiods 
       convert(HarmonicArray,[convert(PeriodicFunctionMatrix,A) convert(PeriodicFunctionMatrix,B)])
    end
end
hcat(A::HarmonicArray, B::HarmonicArray) = horzcat(A,B)
hcat(A::HarmonicArray, C::AbstractMatrix) = horzcat(A, HarmonicArray(C, A.period))
hcat(A::AbstractMatrix, C::HarmonicArray) = horzcat(HarmonicArray(A, C.period), C)
horzcat(A::HarmonicArray, C::AbstractMatrix) = horzcat(A, HarmonicArray(C, A.period))
horzcat(A::AbstractMatrix, C::HarmonicArray) = horzcat(HarmonicArray(A, C.period), C)
horzcat(A::HarmonicArray, C::PeriodicFunctionMatrix) = horzcat(convert(PeriodicFunctionMatrix,A), C)
horzcat(A::PeriodicFunctionMatrix, C::HarmonicArray) = horzcat(A, convert(PeriodicFunctionMatrix,C))
hcat(A::HarmonicArray, C::PeriodicFunctionMatrix) = horzcat(convert(PeriodicFunctionMatrix,A), C)
hcat(A::PeriodicFunctionMatrix, C::HarmonicArray) = horzcat(A, convert(PeriodicFunctionMatrix,C))

function vertcat(A::HarmonicArray, B::HarmonicArray)
    if A.period == B.period && A.nperiod == B.nperiod
       m, n, la = size(A.values)
       mb, nb, lb = size(B.values)
       n == nb || throw(DimensionMismatch("A and B must have the same number of columns"))
       T = promote_type(eltype(A),eltype(B))
       lmax = max(la,lb)
       Ahr = zeros(Complex{T},m+mb,n,lmax)
       copyto!(view(Ahr,1:m,1:n,1:la), view(A.values,1:m,1:n,1:la))
       copyto!(view(Ahr,m+1:m+mb,1:n,1:lb), view(B.values,1:mb,1:n,1:lb))
       return HarmonicArray{:c,T}(Ahr, A.period, nperiod = A.nperiod) 
    else
       #TODO: fix different numbers of subperiods 
       convert(HarmonicArray,[convert(PeriodicFunctionMatrix,A); convert(PeriodicFunctionMatrix,B)])
    end
end
vcat(A::HarmonicArray, B::HarmonicArray) = vertcat(A,B)
vcat(A::HarmonicArray, C::AbstractMatrix) = vertcat(A, HarmonicArray(C, A.period))
vcat(A::AbstractMatrix, C::HarmonicArray) = vertcat(HarmonicArray(A, C.period), C)
vertcat(A::HarmonicArray, C::AbstractMatrix) = vertcat(A, HarmonicArray(C, A.period))
vertcat(A::AbstractMatrix, C::HarmonicArray) = vertcat(HarmonicArray(A, C.period), C)
vertcat(A::HarmonicArray, C::PeriodicFunctionMatrix) = vertcat(convert(PeriodicFunctionMatrix,A), C)
vertcat(A::PeriodicFunctionMatrix, C::HarmonicArray) = vertcat(A, convert(PeriodicFunctionMatrix,C))
vcat(A::HarmonicArray, C::PeriodicFunctionMatrix) = vertcat(convert(PeriodicFunctionMatrix,A), C)
vcat(A::PeriodicFunctionMatrix, C::HarmonicArray) = vertcat(A, convert(PeriodicFunctionMatrix,C))


function blockdiag(A::HarmonicArray, B::HarmonicArray)
    if A.period == B.period && A.nperiod == B.nperiod
        ma, na, la = size(A.values)
        mb, nb, lb = size(B.values)
        T = promote_type(eltype(A),eltype(B))
        lmax = max(la,lb)
        Ahr = zeros(Complex{T},ma+mb,na+nb,lmax)
        copyto!(view(Ahr,1:ma,1:na,1:la),A.values) 
        copyto!(view(Ahr,ma+1:ma+mb,na+1:na+nb,1:lb),B.values) 
        return HarmonicArray{:c,real(T)}(Ahr, A.period, nperiod = A.nperiod) 
    else
        #TODO: fix different numbers of subperiods 
        convert(HarmonicArray,blockdiag(convert(PeriodicFunctionMatrix,A),convert(PeriodicFunctionMatrix,B)))
    end
end


# Operations with periodic time-series matrices
function pmderiv(A::PeriodicTimeSeriesMatrix{:c,T}; kwargs...) where {T}
    N = length(A)
    #tvmdereval(A, (0:N-1)*A.period/A.nperiod/N)
    #PeriodicTimeSeriesMatrix{:c,T}(tvmeval(pmderiv(convert(HarmonicArray,A)), collect((0:N-1)*A.period/A.nperiod/N)), A.period, A.nperiod)
    convert(PeriodicTimeSeriesMatrix,pmderiv(convert(PeriodicFunctionMatrix,A);kwargs...); ns = length(A))
end
LinearAlgebra.inv(A::PeriodicTimeSeriesMatrix) = PeriodicTimeSeriesMatrix(inv.(A.values), A.period; nperiod = A.nperiod)
LinearAlgebra.transpose(A::PeriodicTimeSeriesMatrix{:c,T}) where {T} = 
    PeriodicTimeSeriesMatrix{:c,T}([copy(transpose(A.values[i])) for i in 1:length(A)], A.period, A.nperiod)
LinearAlgebra.adjoint(A::PeriodicTimeSeriesMatrix) = transpose(A)
#LinearAlgebra.tr(A::PeriodicTimeSeriesMatrix) = PeriodicTimeSeriesMatrix(tr.(A.values), A.period; nperiod = A.nperiod)
LinearAlgebra.tr(A::PeriodicTimeSeriesMatrix) = PeriodicTimeSeriesMatrix([[tr(A.values[i])] for i in 1:length(A)], A.period; nperiod = A.nperiod)
function trace(A::PeriodicTimeSeriesMatrix) 
    K = length(A)
    K == 0 && (return zeros(eltype(A)))
    tt = tr(A.values[1])
    for i = 2:K
        tt += tr(A.values[i])
    end 
    return tt*A.nperiod/K
end
LinearAlgebra.eigvals(A::PeriodicTimeSeriesMatrix) = eigvals.(A.values)
function LinearAlgebra.opnorm(A::PeriodicTimeSeriesMatrix, p::Union{Real,Missing} = missing)
    if ismissing(p)
        return PeriodicTimeSeriesMatrix([[norm(A.values[i])] for i in 1:length(A)], A.period; nperiod = A.nperiod)
    else
        return PeriodicTimeSeriesMatrix([[opnorm(A.values[i],p)] for i in 1:length(A)], A.period; nperiod = A.nperiod)
    end
end
function LinearAlgebra.norm(A::PeriodicTimeSeriesMatrix, p::Real = 2) 
    K = length(A)
    K == 0 && (return zeros(eltype(A)))
    nrm = zero(eltype(A))
    if p == 2 
       for i = 1:K
           nrm += norm(A.values[i])^2
       end 
       return sqrt(nrm*A.period/K)
    elseif isinf(p) 
        for i = 1:K
            nrm = max(nrm,norm(A.values[i]))
        end 
        return nrm
    elseif p == 1 
        for i = 1:K
            nrm += norm(A.values[i])
        end 
        return nrm*A.period/K
    else 
       throw(ArgumentError("only p-norms for p = 1, 2, or Inf are supported"))
    end
end
function +(A::PeriodicTimeSeriesMatrix, B::PeriodicTimeSeriesMatrix)
    A.period == B.period && A.nperiod == B.nperiod && length(A) == length(B) &&
        (return PeriodicTimeSeriesMatrix{:c,promote_type(eltype(A),eltype(B))}([A.values[i]+B.values[i] for i in 1:length(A)], A.period, A.nperiod))
    isconstant(A) && 
       (return PeriodicTimeSeriesMatrix{:c,promote_type(eltype(A),eltype(B))}([B.values[i]+A.values[1] for i in 1:length(B)], B.period, B.nperiod))
    isconstant(B) && 
       (return PeriodicTimeSeriesMatrix{:c,promote_type(eltype(A),eltype(B))}([A.values[i]+B.values[1] for i in 1:length(A)], A.period, A.nperiod))
    if A.period == B.period 
       nperiod = gcd(A.nperiod,B.nperiod)
       ns = div(lcm(A.nperiod*length(A),B.nperiod*length(B)),nperiod)
       Δ = A.period/nperiod/ns
       δ = Δ/2
       return PeriodicTimeSeriesMatrix{:c,promote_type(eltype(A),eltype(B))}([tpmeval(A,(i-1)*Δ+δ)+tpmeval(B,(i-1)*Δ+δ) for i in 1:ns], A.period, nperiod) 
    else       
       Tsub = A.period/A.nperiod
       Tsub ≈ B.period/B.nperiod || error("periods or subperiods must be equal for addition")
       nperiod = lcm(A.nperiod,B.nperiod)
       period = Tsub*nperiod
       ns = lcm(length(A),length(B))
       Δ = Tsub/ns
       δ = Δ/2
       return PeriodicTimeSeriesMatrix{:c,promote_type(eltype(A),eltype(B))}([tpmeval(A,(i-1)*Δ+δ)+tpmeval(B,(i-1)*Δ+δ) for i in 1:ns], period, nperiod)   
    end     
end
+(A::PeriodicTimeSeriesMatrix, C::AbstractMatrix) = +(A, PeriodicTimeSeriesMatrix(C, A.period))
+(A::AbstractMatrix, C::PeriodicTimeSeriesMatrix) = +(PeriodicTimeSeriesMatrix(A, C.period), C)
+(A::PeriodicTimeSeriesMatrix, C::PeriodicFunctionMatrix) = +(convert(PeriodicFunctionMatrix,A;method="constant"), C)
+(A::PeriodicFunctionMatrix, C::PeriodicTimeSeriesMatrix) = +(A, convert(PeriodicFunctionMatrix,C;method="constant"))
-(A::PeriodicTimeSeriesMatrix) = PeriodicTimeSeriesMatrix{:c,eltype(A)}([-A.values[i] for i in 1:length(A)], A.period, A.nperiod)
-(A::PeriodicTimeSeriesMatrix, B::PeriodicTimeSeriesMatrix) = +(A,-B)
-(A::PeriodicTimeSeriesMatrix, C::AbstractMatrix) = +(A,-C)
-(A::AbstractMatrix, C::PeriodicTimeSeriesMatrix) = +(A, -C)
-(A::PeriodicTimeSeriesMatrix, C::PeriodicFunctionMatrix) = +(A,-C)
-(A::PeriodicFunctionMatrix, C::PeriodicTimeSeriesMatrix) = +(A, -C)

function (+)(A::PeriodicTimeSeriesMatrix, J::UniformScaling{<:Real}) 
    m, n = size(A)
    n == m || throw(DimensionMismatch("matrix is not square: dimensions are $((m,n))"))
    A+Matrix(J(n))
end
(+)(J::UniformScaling{<:Real}, A::PeriodicTimeSeriesMatrix) = +(A,J)
(-)(A::PeriodicTimeSeriesMatrix, J::UniformScaling{<:Real}) = +(A,-J)
(-)(J::UniformScaling{<:Real}, A::PeriodicTimeSeriesMatrix) = +(-A,J)

function *(A::PeriodicTimeSeriesMatrix, B::PeriodicTimeSeriesMatrix)
    A.period == B.period && A.nperiod == B.nperiod && length(A) == length(B) &&
        (return PeriodicTimeSeriesMatrix{:c,promote_type(eltype(A),eltype(B))}([A.values[i]*B.values[i] for i in 1:length(A)], A.period, A.nperiod))
    isconstant(A) && 
       (return PeriodicTimeSeriesMatrix{:c,promote_type(eltype(A),eltype(B))}([A.values[1]*B.values[i] for i in 1:length(B)], B.period, B.nperiod))
    isconstant(B) && 
       (return PeriodicTimeSeriesMatrix{:c,promote_type(eltype(A),eltype(B))}([A.values[i]*B.values[1] for i in 1:length(A)], A.period, A.nperiod))
    if A.period == B.period 
       nperiod = gcd(A.nperiod,B.nperiod)
       ns = div(lcm(A.nperiod*length(A),B.nperiod*length(B)),nperiod)
       Δ = A.period/nperiod/ns
       δ = Δ/2
       return PeriodicTimeSeriesMatrix{:c,promote_type(eltype(A),eltype(B))}([tpmeval(A,(i-1)*Δ+δ)*tpmeval(B,(i-1)*Δ+δ) for i in 1:ns], A.period, nperiod) 
    else          
       Tsub = A.period/A.nperiod
       Tsub ≈ B.period/B.nperiod || error("periods or subperiods must be equal for multiplication")
       nperiod = lcm(A.nperiod,B.nperiod)
       period = Tsub*nperiod
       ns = lcm(length(A),length(B))
       Δ = Tsub/ns
       δ = Δ/2
       return PeriodicTimeSeriesMatrix{:c,promote_type(eltype(A),eltype(B))}([tpmeval(A,(i-1)*Δ+δ)*tpmeval(B,(i-1)*Δ+δ) for i in 1:ns], period, nperiod)   
    end     
end
*(A::PeriodicTimeSeriesMatrix, C::AbstractMatrix) = *(A, PeriodicTimeSeriesMatrix(C, A.period))
*(A::AbstractMatrix, C::PeriodicTimeSeriesMatrix) = *(PeriodicTimeSeriesMatrix(A, C.period), C)
*(A::PeriodicTimeSeriesMatrix, C::Real) = PeriodicTimeSeriesMatrix{:c,promote_type(eltype(A),eltype(C))}([A.values[i]*C for i in 1:length(A)], A.period, A.nperiod)
*(C::Real, A::PeriodicTimeSeriesMatrix) = PeriodicTimeSeriesMatrix{:c,promote_type(eltype(A),eltype(C))}([C*A.values[i] for i in 1:length(A)], A.period, A.nperiod)
*(A::PeriodicTimeSeriesMatrix, C::PeriodicFunctionMatrix) = *(convert(PeriodicFunctionMatrix,A;method="constant"), C)
*(A::PeriodicFunctionMatrix, C::PeriodicTimeSeriesMatrix) = *(A, convert(PeriodicFunctionMatrix,C;method="constant"))
/(A::PeriodicTimeSeriesMatrix, C::Real) = *(A, 1/C)
*(J::UniformScaling{<:Real}, A::PeriodicTimeSeriesMatrix) = J.λ*A
*(A::PeriodicTimeSeriesMatrix, J::UniformScaling{<:Real}) = A*J.λ


Base.iszero(A::PeriodicTimeSeriesMatrix) = iszero(A.values)
LinearAlgebra.issymmetric(A::PeriodicTimeSeriesMatrix) = all(issymmetric.(A.values))
function ==(A::PeriodicTimeSeriesMatrix, B::PeriodicTimeSeriesMatrix)
    A.period == B.period && A.nperiod == B.nperiod && isequal(A.values, B.values)
end
function Base.isapprox(A::PeriodicTimeSeriesMatrix, B::PeriodicTimeSeriesMatrix; rtol::Real = sqrt(eps(Float64)), atol::Real = 0)
    A.period == B.period && A.nperiod == B.nperiod && isapprox(A.values, B.values; rtol, atol)
end
function Base.isapprox(A::PeriodicTimeSeriesMatrix, J::UniformScaling{<:Real}; kwargs...)
    all([isapprox(A.values[i], J; kwargs...) for i in 1:length(A.values)])
end
Base.isapprox(J::UniformScaling{<:Real}, A::PeriodicTimeSeriesMatrix; kwargs...) = isapprox(A, J; kwargs...)

function horzcat(A::PeriodicTimeSeriesMatrix, B::PeriodicTimeSeriesMatrix)
    A.period == B.period && A.nperiod == B.nperiod && length(A) == length(B) && 
        (return PeriodicTimeSeriesMatrix{:c,promote_type(eltype(A),eltype(B))}([[A.values[i] B.values[i]] for i in 1:length(A)], A.period, A.nperiod))
    isconstant(A) && 
       (return PeriodicTimeSeriesMatrix{:c,promote_type(eltype(A),eltype(B))}([[A.values[1] B.values[i]] for i in 1:length(B)], B.period, B.nperiod))
    isconstant(B) && 
       (return PeriodicTimeSeriesMatrix{:c,promote_type(eltype(A),eltype(B))}([[A.values[i] B.values[1]] for i in 1:length(A)], A.period, A.nperiod))
    if A.period == B.period 
        nperiod = gcd(A.nperiod,B.nperiod)
        ns = div(lcm(A.nperiod*length(A),B.nperiod*length(B)),nperiod)
        Δ = A.period/nperiod/ns
        δ = Δ/2
        return PeriodicTimeSeriesMatrix{:c,promote_type(eltype(A),eltype(B))}([[tpmeval(A,(i-1)*Δ+δ) tpmeval(B,(i-1)*Δ+δ)] for i in 1:ns], A.period, nperiod) 
    else          
        Tsub = A.period/A.nperiod
        Tsub ≈ B.period/B.nperiod || error("periods or subperiods must be equal for horizontal concatenation")
        nperiod = lcm(A.nperiod,B.nperiod)
        period = Tsub*nperiod
        ns = lcm(length(A),length(B))
        Δ = Tsub/ns
        δ = Δ/2
        return PeriodicTimeSeriesMatrix{:c,promote_type(eltype(A),eltype(B))}([[tpmeval(A,(i-1)*Δ+δ) tpmeval(B,(i-1)*Δ+δ)] for i in 1:ns], period, nperiod)   
    end     
end
hcat(A::PeriodicTimeSeriesMatrix, B::PeriodicTimeSeriesMatrix) = horzcat(A,B)
hcat(A::PeriodicTimeSeriesMatrix, C::AbstractMatrix) = horzcat(A, PeriodicTimeSeriesMatrix(C, A.period))
hcat(A::AbstractMatrix, C::PeriodicTimeSeriesMatrix) = horzcat(PeriodicTimeSeriesMatrix(A, C.period), C)
horzcat(A::PeriodicTimeSeriesMatrix, C::AbstractMatrix) = horzcat(A, PeriodicTimeSeriesMatrix(C, A.period))
horzcat(A::AbstractMatrix, C::PeriodicTimeSeriesMatrix) = horzcat(PeriodicTimeSeriesMatrix(A, C.period), C)


function vertcat(A::PeriodicTimeSeriesMatrix, B::PeriodicTimeSeriesMatrix)
    A.period == B.period && A.nperiod == B.nperiod && length(A) == length(B) && 
        (return PeriodicTimeSeriesMatrix{:c,promote_type(eltype(A),eltype(B))}([[A.values[i]; B.values[i]] for i in 1:length(A)], A.period, A.nperiod))
    isconstant(A) && 
       (return PeriodicTimeSeriesMatrix{:c,promote_type(eltype(A),eltype(B))}([[A.values[1]; B.values[i]] for i in 1:length(B)], B.period, B.nperiod))
    isconstant(B) && 
       (return PeriodicTimeSeriesMatrix{:c,promote_type(eltype(A),eltype(B))}([[A.values[i]; B.values[1]] for i in 1:length(A)], A.period, A.nperiod))
    if A.period == B.period 
        nperiod = gcd(A.nperiod,B.nperiod)
        ns = div(lcm(A.nperiod*length(A),B.nperiod*length(B)),nperiod)
        Δ = A.period/nperiod/ns
        δ = Δ/2
        return PeriodicTimeSeriesMatrix{:c,promote_type(eltype(A),eltype(B))}([[tpmeval(A,(i-1)*Δ+δ); tpmeval(B,(i-1)*Δ+δ)] for i in 1:ns], A.period, nperiod) 
    else          
        Tsub = A.period/A.nperiod
        Tsub ≈ B.period/B.nperiod || error("periods or subperiods must be equal for vertical concatenation")
        nperiod = lcm(A.nperiod,B.nperiod)
        period = Tsub*nperiod
        ns = lcm(length(A),length(B))
        Δ = Tsub/ns
        δ = Δ/2
        return PeriodicTimeSeriesMatrix{:c,promote_type(eltype(A),eltype(B))}([[tpmeval(A,(i-1)*Δ+δ); tpmeval(B,(i-1)*Δ+δ)] for i in 1:ns], period, nperiod)   
    end     
end
vcat(A::PeriodicTimeSeriesMatrix, B::PeriodicTimeSeriesMatrix) = vertcat(A,B)
vcat(A::PeriodicTimeSeriesMatrix, C::AbstractMatrix) = vertcat(A, PeriodicTimeSeriesMatrix(C, A.period))
vcat(A::AbstractMatrix, C::PeriodicTimeSeriesMatrix) = vertcat(PeriodicTimeSeriesMatrix(A, C.period), C)
vertcat(A::PeriodicTimeSeriesMatrix, C::AbstractMatrix) = vertcat(A, PeriodicTimeSeriesMatrix(C, A.period))
vertcat(A::AbstractMatrix, C::PeriodicTimeSeriesMatrix) = vertcat(PeriodicTimeSeriesMatrix(A, C.period), C)

function blockdiag(A::PeriodicTimeSeriesMatrix, B::PeriodicTimeSeriesMatrix)
    A.period == B.period && A.nperiod == B.nperiod && length(A) == length(B) && 
        (return PeriodicTimeSeriesMatrix{:c,promote_type(eltype(A),eltype(B))}([bldiag(A.values[i],B.values[i]) for i in 1:length(A)], A.period, A.nperiod))
    isconstant(A) && 
       (return PeriodicTimeSeriesMatrix{:c,promote_type(eltype(A),eltype(B))}([bldiag(A.values[1],B.values[i]) for i in 1:length(B)], B.period, B.nperiod))
    isconstant(B) && 
       (return PeriodicTimeSeriesMatrix{:c,promote_type(eltype(A),eltype(B))}([bldiag(A.values[i],B.values[1]) for i in 1:length(A)], A.period, A.nperiod))
    if A.period == B.period 
        nperiod = gcd(A.nperiod,B.nperiod)
        ns = div(lcm(A.nperiod*length(A),B.nperiod*length(B)),nperiod)
        Δ = A.period/nperiod/ns
        δ = Δ/2
        return PeriodicTimeSeriesMatrix{:c,promote_type(eltype(A),eltype(B))}([bldiag(tpmeval(A,(i-1)*Δ+δ), tpmeval(B,(i-1)*Δ+δ)) for i in 1:ns], A.period, nperiod) 
    else          
        Tsub = A.period/A.nperiod
        Tsub ≈ B.period/B.nperiod || error("periods or subperiods must be equal for block-diagonal appending")
        nperiod = lcm(A.nperiod,B.nperiod)
        period = Tsub*nperiod
        ns = lcm(length(A),length(B))
        Δ = Tsub/ns
        δ = Δ/2
        return PeriodicTimeSeriesMatrix{:c,promote_type(eltype(A),eltype(B))}([bldiag(tpmeval(A,(i-1)*Δ+δ),tpmeval(B,(i-1)*Δ+δ)) for i in 1:ns], period, nperiod)   
    end     
end



# Operations with periodic switching matrices
# function pmderiv(A::PeriodicSwitchingMatrix{:c,T}) where {T}
#     PeriodicSwitchingMatrix{:c,T}([zeros(T,size(A,1),size(A,2)) for i in 1:length(A)], A.ts, A.period, A.nperiod)
# end
function pmderiv(A::PeriodicSwitchingMatrix{:c,T}; kwargs...) where {T}
    @warn "No derivative function available for a periodic switching matrix"
    #convert(PeriodicSwitchingMatrix,pmderiv(convert(PeriodicFunctionMatrix,A);kwargs...))
end

LinearAlgebra.inv(A::PeriodicSwitchingMatrix) = PeriodicSwitchingMatrix(inv.(A.values), A.ts, A.period; nperiod = A.nperiod)
LinearAlgebra.transpose(A::PeriodicSwitchingMatrix{:c,T}) where {T} = 
    PeriodicSwitchingMatrix{:c,T}([copy(transpose(A.values[i])) for i in 1:length(A)], A.ts, A.period, A.nperiod)
LinearAlgebra.adjoint(A::PeriodicSwitchingMatrix) = transpose(A)
#LinearAlgebra.tr(A::PeriodicSwitchingMatrix) = PeriodicSwitchingMatrix(tr.(A.values), A.period; nperiod = A.nperiod)
LinearAlgebra.tr(A::PeriodicSwitchingMatrix) = PeriodicSwitchingMatrix([[tr(A.values[i])] for i in 1:length(A)], A.ts, A.period; nperiod = A.nperiod)
function trace(A::PeriodicSwitchingMatrix) 
    K = length(A)
    K == 0 && (return zeros(eltype(A)))
    tt = tr(A.values[K])*(A.period/A.nperiod-A.ts[K])
    for i = 1:K-1
        tt += tr(A.values[i])*(A.ts[i+1]-A.ts[i])
    end 
    return tt*A.nperiod/A.period
end

LinearAlgebra.eigvals(A::PeriodicSwitchingMatrix) = eigvals.(A.values)
function LinearAlgebra.opnorm(A::PeriodicSwitchingMatrix, p::Union{Real,Missing} = missing)
    if ismissing(p)
        return PeriodicSwitchingMatrix([[norm(A.values[i])] for i in 1:length(A)], A.ts, A.period; nperiod = A.nperiod)
    else
        return PeriodicSwitchingMatrix([[opnorm(A.values[i],p)] for i in 1:length(A)], A.ts, A.period; nperiod = A.nperiod)
    end
end
function LinearAlgebra.norm(A::PeriodicSwitchingMatrix, p::Real = 2) 
    K = length(A)
    K == 0 && (return zeros(eltype(A)))
    Δ = A.period/A.nperiod/K
    ts = zero(eltype(Δ))
    if p == 2
        nrm = norm(A.values[K])^2*(A.period/A.nperiod-A.ts[K])
        for i = 1:K-1
            nrm += norm(A.values[i])^2*(A.ts[i+1]-A.ts[i])
        end 
        return sqrt(nrm*A.nperiod)
    elseif isinf(p)
        nrm = zero(eltype(A))
        for i = 1:K-1
            nrm = max(nrm,norm(A.values[i]))
        end 
        return nrm
    elseif p == 1    
        nrm = norm(A.values[K])*(A.period/A.nperiod-A.ts[K])
        for i = 1:K-1
            nrm += norm(A.values[i])*(A.ts[i+1]-A.ts[i])
        end 
        return nrm*A.nperiod
    else
        throw(ArgumentError("only p-norms for p = 1, 2, or Inf are supported"))
    end
end
function +(A::PeriodicSwitchingMatrix, B::PeriodicSwitchingMatrix)
    A.period == B.period && A.nperiod == B.nperiod && length(A) == length(B) && A.ts == B.ts &&
        (return PeriodicSwitchingMatrix{:c,promote_type(eltype(A),eltype(B))}([A.values[i]+B.values[i] for i in 1:length(A)], A.ts, A.period, A.nperiod))
    isconstant(A) && 
       (return PeriodicSwitchingMatrix{:c,promote_type(eltype(A),eltype(B))}([B.values[i]+A.values[1] for i in 1:length(B)], B.ts, B.period, B.nperiod))
    isconstant(B) && 
       (return PeriodicSwitchingMatrix{:c,promote_type(eltype(A),eltype(B))}([A.values[i]+B.values[1] for i in 1:length(A)], A.ts, A.period, A.nperiod))
    A.period == B.period || error("periods must be equal for addition")
    if A.nperiod == B.nperiod
        ts = unique(sort([A.ts;B.ts]))
    else
        ts = unique(sort([vcat([(i-1)*A.period/A.nperiod .+ A.ts for i in 1:A.nperiod]...);
                          vcat([(i-1)*B.period/B.nperiod .+ B.ts for i in 1:B.nperiod]...)]))
    end
    return PeriodicSwitchingMatrix{:c,promote_type(eltype(A),eltype(B))}([tpmeval(A,ts[i])+tpmeval(B,ts[i]) for i in 1:length(ts)], ts, A.period, gcd(A.nperiod,B.nperiod))
end
+(A::PeriodicSwitchingMatrix, C::AbstractMatrix) = +(A, PeriodicSwitchingMatrix([C], [0], A.period))
+(A::AbstractMatrix, C::PeriodicSwitchingMatrix) = +(PeriodicSwitchingMatrix([A], [0], C.period), C)
+(A::PeriodicSwitchingMatrix, C::PeriodicFunctionMatrix) = +(convert(PeriodicFunctionMatrix,A), C)
+(A::PeriodicFunctionMatrix, C::PeriodicSwitchingMatrix) = +(A, convert(PeriodicFunctionMatrix,C))
-(A::PeriodicSwitchingMatrix) = PeriodicSwitchingMatrix{:c,eltype(A)}([-A.values[i] for i in 1:length(A)], A.ts, A.period, A.nperiod)
-(A::PeriodicSwitchingMatrix, B::PeriodicSwitchingMatrix) = +(A,-B)
-(A::PeriodicSwitchingMatrix, C::AbstractMatrix) = +(A,-C)
-(A::AbstractMatrix, C::PeriodicSwitchingMatrix) = +(A, -C)
-(A::PeriodicSwitchingMatrix, C::PeriodicFunctionMatrix) = +(A,-C)
-(A::PeriodicFunctionMatrix, C::PeriodicSwitchingMatrix) = +(A, -C)

function (+)(A::PeriodicSwitchingMatrix, J::UniformScaling{<:Real}) 
    m, n = size(A)
    n == m || throw(DimensionMismatch("matrix is not square: dimensions are $((m,n))"))
    A+Matrix(J(n))
end
(+)(J::UniformScaling{<:Real}, A::PeriodicSwitchingMatrix) = +(A,J)
(-)(A::PeriodicSwitchingMatrix, J::UniformScaling{<:Real}) = +(A,-J)
(-)(J::UniformScaling{<:Real}, A::PeriodicSwitchingMatrix) = +(-A,J)

function *(A::PeriodicSwitchingMatrix, B::PeriodicSwitchingMatrix)
    A.period == B.period && A.nperiod == B.nperiod && length(A) == length(B) && A.ts == B.ts &&
        (return PeriodicSwitchingMatrix{:c,promote_type(eltype(A),eltype(B))}([A.values[i]*B.values[i] for i in 1:length(A)], A.ts, A.period, A.nperiod))
    isconstant(A) && 
       (return PeriodicSwitchingMatrix{:c,promote_type(eltype(A),eltype(B))}([A.values[1]*B.values[i] for i in 1:length(B)], B.ts, B.period, B.nperiod))
    isconstant(B) && 
       (return PeriodicSwitchingMatrix{:c,promote_type(eltype(A),eltype(B))}([A.values[i]*B.values[1] for i in 1:length(A)], A.ts, A.period, A.nperiod))
    A.period == B.period || error("periods must be equal for multiplication")
    if A.nperiod == B.nperiod
        ts = unique(sort([A.ts;B.ts]))
    else
        ts = unique(sort([vcat([(i-1)*A.period/A.nperiod .+ A.ts for i in 1:A.nperiod]...);
                          vcat([(i-1)*B.period/B.nperiod .+ B.ts for i in 1:B.nperiod]...)]))
    end
    return PeriodicSwitchingMatrix{:c,promote_type(eltype(A),eltype(B))}([tpmeval(A,ts[i])*tpmeval(B,ts[i]) for i in 1:length(ts)], ts, A.period, gcd(A.nperiod,B.nperiod))
   end
*(A::PeriodicSwitchingMatrix, C::AbstractMatrix) = PeriodicSwitchingMatrix{:c,promote_type(eltype(A),eltype(C))}([A.values[i]*C for i in 1:length(A)], A.ts, A.period, A.nperiod)
*(A::AbstractMatrix, C::PeriodicSwitchingMatrix) = PeriodicSwitchingMatrix{:c,promote_type(eltype(A),eltype(C))}([A*C.values[i] for i in 1:length(C)], C.ts, C.period, C.nperiod)
*(A::PeriodicSwitchingMatrix, C::Real) = PeriodicSwitchingMatrix{:c,promote_type(eltype(A),eltype(C))}([A.values[i]*C for i in 1:length(A)], A.ts, A.period, A.nperiod)
*(C::Real, A::PeriodicSwitchingMatrix) = PeriodicSwitchingMatrix{:c,promote_type(eltype(A),eltype(C))}([C*A.values[i] for i in 1:length(A)], A.ts, A.period, A.nperiod)
*(A::PeriodicSwitchingMatrix, C::PeriodicFunctionMatrix) = *(convert(PeriodicFunctionMatrix,A), C)
*(A::PeriodicFunctionMatrix, C::PeriodicSwitchingMatrix) = *(A, convert(PeriodicFunctionMatrix,C))
/(A::PeriodicSwitchingMatrix, C::Real) = *(A, 1/C)
*(J::UniformScaling{<:Real}, A::PeriodicSwitchingMatrix) = J.λ*A
*(A::PeriodicSwitchingMatrix, J::UniformScaling{<:Real}) = A*J.λ

Base.iszero(A::PeriodicSwitchingMatrix) = iszero(A.values)
LinearAlgebra.issymmetric(A::PeriodicSwitchingMatrix) = all(issymmetric.(A.values))
function ==(A::PeriodicSwitchingMatrix, B::PeriodicSwitchingMatrix)
    A.period == B.period && A.nperiod == B.nperiod && A.ts == B.ts && isequal(A.values, B.values)
end
function Base.isapprox(A::PeriodicSwitchingMatrix, B::PeriodicSwitchingMatrix; rtol::Real = sqrt(eps(Float64)), atol::Real = 0)
    A.period == B.period && A.nperiod == B.nperiod && A.ts == B.ts && isapprox(A.values, B.values; rtol, atol) 
end
function Base.isapprox(A::PeriodicSwitchingMatrix, J::UniformScaling{<:Real}; kwargs...)
    all([isapprox(A.values[i], J; kwargs...) for i in 1:length(A.values)])
end
Base.isapprox(J::UniformScaling{<:Real}, A::PeriodicSwitchingMatrix; kwargs...) = isapprox(A, J; kwargs...)

function horzcat(A::PeriodicSwitchingMatrix, B::PeriodicSwitchingMatrix)
    A.period == B.period && A.nperiod == B.nperiod && length(A) == length(B) && A.ts == B.ts &&
        (return PeriodicSwitchingMatrix{:c,promote_type(eltype(A),eltype(B))}([[A.values[i] B.values[i]] for i in 1:length(A)], A.ts, A.period, A.nperiod))
    isconstant(A) && 
       (return PeriodicSwitchingMatrix{:c,promote_type(eltype(A),eltype(B))}([[A.values[1] B.values[i]] for i in 1:length(B)], B.ts, B.period, B.nperiod))
    isconstant(B) && 
       (return PeriodicSwitchingMatrix{:c,promote_type(eltype(A),eltype(B))}([[A.values[i] B.values[1]] for i in 1:length(A)], A.ts, A.period, A.nperiod))
    A.period == B.period || error("periods must be equal for horizontal concatenation")
    if A.nperiod == B.nperiod
        ts = unique(sort([A.ts;B.ts]))
    else
        ts = unique(sort([vcat([(i-1)*A.period/A.nperiod .+ A.ts for i in 1:A.nperiod]...);
                          vcat([(i-1)*B.period/B.nperiod .+ B.ts for i in 1:B.nperiod]...)]))
    end
    return PeriodicSwitchingMatrix{:c,promote_type(eltype(A),eltype(B))}([[tpmeval(A,ts[i]) tpmeval(B,ts[i])] for i in 1:length(ts)], ts, A.period, gcd(A.nperiod,B.nperiod))
end
hcat(A::PeriodicSwitchingMatrix, B::PeriodicSwitchingMatrix) = horzcat(A,B)
hcat(A::PeriodicSwitchingMatrix, C::AbstractMatrix) = horzcat(A, PeriodicSwitchingMatrix(C, A.period))
hcat(A::AbstractMatrix, C::PeriodicSwitchingMatrix) = horzcat(PeriodicSwitchingMatrix(A, C.period), C)
horzcat(A::PeriodicSwitchingMatrix, C::AbstractMatrix) = horzcat(A, PeriodicSwitchingMatrix(C, A.period))
horzcat(A::AbstractMatrix, C::PeriodicSwitchingMatrix) = horzcat(PeriodicSwitchingMatrix(A, C.period), C)

function vertcat(A::PeriodicSwitchingMatrix, B::PeriodicSwitchingMatrix)
    A.period == B.period && A.nperiod == B.nperiod && length(A) == length(B) && A.ts == B.ts &&
        (return PeriodicSwitchingMatrix{:c,promote_type(eltype(A),eltype(B))}([[A.values[i]; B.values[i]] for i in 1:length(A)], A.ts, A.period, A.nperiod))
    isconstant(A) && 
       (return PeriodicSwitchingMatrix{:c,promote_type(eltype(A),eltype(B))}([[A.values[1]; B.values[i]] for i in 1:length(B)], B.ts, B.period, B.nperiod))
    isconstant(B) && 
       (return PeriodicSwitchingMatrix{:c,promote_type(eltype(A),eltype(B))}([[A.values[i]; B.values[1]] for i in 1:length(A)], A.ts, A.period, A.nperiod))
    A.period == B.period || error("periods must be equal for vertical concatenation")
    if A.nperiod == B.nperiod
        ts = unique(sort([A.ts;B.ts]))
    else
        ts = unique(sort([vcat([(i-1)*A.period/A.nperiod .+ A.ts for i in 1:A.nperiod]...);
                          vcat([(i-1)*B.period/B.nperiod .+ B.ts for i in 1:B.nperiod]...)]))
    end
    return PeriodicSwitchingMatrix{:c,promote_type(eltype(A),eltype(B))}([[tpmeval(A,ts[i]); tpmeval(B,ts[i])] for i in 1:length(ts)], ts, A.period, gcd(A.nperiod,B.nperiod))
end
vcat(A::PeriodicSwitchingMatrix, B::PeriodicSwitchingMatrix) = vertcat(A,B)
vcat(A::PeriodicSwitchingMatrix, C::AbstractMatrix) = vertcat(A, PeriodicSwitchingMatrix(C, A.period))
vcat(A::AbstractMatrix, C::PeriodicSwitchingMatrix) = vertcat(PeriodicSwitchingMatrix(A, C.period), C)
vertcat(A::PeriodicSwitchingMatrix, C::AbstractMatrix) = vertcat(A, PeriodicSwitchingMatrix(C, A.period))
vertcat(A::AbstractMatrix, C::PeriodicSwitchingMatrix) = vertcat(PeriodicSwitchingMatrix(A, C.period), C)


function blockdiag(A::PeriodicSwitchingMatrix, B::PeriodicSwitchingMatrix)
    A.period == B.period && A.nperiod == B.nperiod && length(A) == length(B) && A.ts == B.ts &&
        (return PeriodicSwitchingMatrix{:c,promote_type(eltype(A),eltype(B))}([bldiag(A.values[i],B.values[i]) for i in 1:length(A)], A.ts, A.period, A.nperiod))
    isconstant(A) && 
       (return PeriodicSwitchingMatrix{:c,promote_type(eltype(A),eltype(B))}([bldiag(A.values[1], B.values[i]) for i in 1:length(B)], B.ts, B.period, B.nperiod))
    isconstant(B) && 
       (return PeriodicSwitchingMatrix{:c,promote_type(eltype(A),eltype(B))}([bldiag(A.values[i], B.values[1]) for i in 1:length(A)], A.ts, A.period, A.nperiod))
    A.period == B.period || error("periods must be equal for block-diagonal appending")
    if A.nperiod == B.nperiod
        ts = unique(sort([A.ts;B.ts]))
    else
        ts = unique(sort([vcat([(i-1)*A.period/A.nperiod .+ A.ts for i in 1:A.nperiod]...);
                          vcat([(i-1)*B.period/B.nperiod .+ B.ts for i in 1:B.nperiod]...)]))
    end
    return PeriodicSwitchingMatrix{:c,promote_type(eltype(A),eltype(B))}([bldiag(tpmeval(A,ts[i]), tpmeval(B,ts[i])) for i in 1:length(ts)], ts, A.period, gcd(A.nperiod,B.nperiod))
end
"""
    blockut(A, B, C) -> D

Block upper triangular appending of periodic matrices `A`, `B` and `C`. 
The resulting `D` is constructed as `D = [A B; 0 C]`. 
`A`, `B` and `C` may have different, but commensurate periods.
"""
function blockut(A11::PeriodicMatrix, A12::PeriodicMatrix, A22::PeriodicMatrix)
    T = promote_type(eltype(A11),eltype(A12),eltype(A22))
    period = promote_period(A11, A12, A22)
    pa11 = length(A11.M)
    pa12 = length(A12.M)
    pa22 = length(A22.M)
    ma11, na11 = size(A11)
    ma12, na12 = size(A12)
    ma22, na22 = size(A22)
    p = lcm(pa11,pa12,pa22)
    nta11 = numerator(rationalize(period/A11.period))
    nta12 = numerator(rationalize(period/A12.period))
    nta22 = numerator(rationalize(period/A22.period))
    K11 = nta11*A11.nperiod*pa11
    K12 = nta12*A12.nperiod*pa12
    K22 = nta22*A22.nperiod*pa22
    K = max(K11,K12,K22)
    X = Vector{Matrix{T}}(undef, p)
    for i = 1:p
        ia11 = mod(i-1,pa11)+1
        ia12 = mod(i-1,pa12)+1
        ia22 = mod(i-1,pa22)+1
        ma11[ia11] == ma12[ia12] || throw(DimensionMismatch("A11 and A12 must have the same number of rows")) 
        na12[ia12] == na22[ia22] || throw(DimensionMismatch("A12 and A22 must have the same number of columns")) 
        X[i] = [A11.M[ia11] A12.M[ia12]; zeros(T,ma22[ia22],na11[ia11]) A22.M[ia22]]
    end
    return PeriodicMatrix(X, period; nperiod = div(K,p))
 end
 function blockut(A11::PeriodicArray, A12::PeriodicArray, A22::PeriodicArray)
    T = promote_type(eltype(A11),eltype(A12),eltype(A22))
    period = promote_period(A11, A12, A22)
    ma11, na11, pa11 = size(A11.M)
    ma12, na12, pa12 = size(A12.M)
    ma11 == ma12 || throw(DimensionMismatch("A11 and A12 must have the same number of rows")) 
    ma22, na22, pa22 = size(A22.M)
    na12 == na22 || throw(DimensionMismatch("A12 and A22 must have the same number of columns")) 
    p = lcm(pa11,pa12,pa22)
    nta11 = numerator(rationalize(period/A11.period))
    nta12 = numerator(rationalize(period/A12.period))
    nta22 = numerator(rationalize(period/A22.period))
    K11 = nta11*A11.nperiod*pa11
    K12 = nta12*A12.nperiod*pa12
    K22 = nta22*A22.nperiod*pa22
    K = max(K11,K12,K22)
    X = Array{T,3}(undef, ma11+ma22, na11+na12, p)
    for i = 1:p
        ia11 = mod(i-1,pa11)+1
        ia12 = mod(i-1,pa12)+1
        ia22 = mod(i-1,pa22)+1
        X[:,:,i] = [[view(A11.M,:,:,ia11) view(A12.M,:,:,ia12)]; [zeros(T,ma22,na11) view(A22.M,:,:,ia22)]]
    end
    return PeriodicArray(X, period; nperiod = div(K,p))
 end
 function blockut(A11::HarmonicArray, A12::HarmonicArray, A22::HarmonicArray)
    if (A11.period == A12.period == A22.period) && (A11.nperiod == A12.nperiod == A22.nperiod)
       ma11, na11, la11 = size(A11.values)
       ma12, na12, la12 = size(A12.values)
       ma22, na22, la22 = size(A22.values) 
       ma11 == ma12 || throw(DimensionMismatch("A11 and A12 must have the same number of rows")) 
       na12 == na22 || throw(DimensionMismatch("A12 and A22 must have the same number of columns"))   
       T = promote_type(eltype(A11),eltype(A12),eltype(A22))
       lmax = max(la11,la12,la22)
       Ahr = zeros(Complex{T},ma11+ma22,na11+na12,lmax)
       copyto!(view(Ahr,1:ma11,1:na11,1:la11),A11.values) 
       copyto!(view(Ahr,1:ma11,na11+1:na11+na12,1:la12),A12.values) 
       copyto!(view(Ahr,ma11+1:ma11+ma22,na11+1:na11+na12,1:la22),A22.values) 
       return HarmonicArray{:c,real(T)}(Ahr, A11.period, nperiod = A11.nperiod) 
    else
       convert(HarmonicArray,blockut(convert(PeriodicFunctionMatrix,A11),convert(PeriodicFunctionMatrix,A12),convert(PeriodicFunctionMatrix,A22)))
    end
 end
 function blockut(A11::PeriodicFunctionMatrix, A12::PeriodicFunctionMatrix, A22::PeriodicFunctionMatrix)
    ma11, na11 = size(A11)
    ma12, na12 = size(A12)
    ma22, na22 = size(A22)   
    ma11 == ma12 || throw(DimensionMismatch("A11 and A12 must have the same number of rows")) 
    na12 == na22 || throw(DimensionMismatch("A12 and A22 must have the same number of columns"))   
    T = promote_type(eltype(A11),eltype(A12),eltype(A22))
    # period = promote_period(A11, A12, A22)
    # nperiod = gcd(A11.nperiod,A12.nperiod,A22.nperiod)
    period, nperiod = promote_period2(A11, A12, A22)
    if isconstant(A11) && isconstant(A12) && isconstant(A22)
       return PeriodicFunctionMatrix{:c,T}(t -> [A11.f(0) A12.f(0); zeros(T,ma22,na11) A22.f(0)] , period, (ma11+ma22,na11+na22), nperiod, true)
    else
       return PeriodicFunctionMatrix{:c,T}(t -> [tpmeval(A11,t) tpmeval(A12,t); zeros(T,ma22,na11) tpmeval(A22,t)], period, (ma11+ma22,na11+na22), nperiod, false)
    end
 end
 
