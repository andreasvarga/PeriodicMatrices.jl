#FourierFunctionMatrices
pmderiv(A::FourierFunctionMatrix{:c,T}) where {T} = FourierFunctionMatrix{:c,T,Fun}(differentiate(A.M), A.period, A.nperiod)
LinearAlgebra.inv(A::FourierFunctionMatrix) = FourierFunctionMatrix(inv(A.M), A.period)
LinearAlgebra.transpose(A::FourierFunctionMatrix{:c,T}) where {T}  = FourierFunctionMatrix{:c,T,Fun}(transpose(A.M), A.period, A.nperiod)
LinearAlgebra.adjoint(A::FourierFunctionMatrix) = transpose(A)
function LinearAlgebra.tr(V::Fun)
    typeof(size(space(V))) == Tuple{} && (return V)
    n, m = size(space(V))
    if n ≠ m
        throw(DimensionMismatch("space $(space(V)) is not square"))
    end
    a = Array(V)
    temp = a[1,1]
    for i = 2:n
        temp += a[i,i]
    end
    return temp
end
LinearAlgebra.tr(A::FourierFunctionMatrix) = FourierFunctionMatrix([tr(A.M);], A.period)
function trace(A::FourierFunctionMatrix; rtol = sqrt(eps())) 
    isconstant(A) && (return tr(tpmeval(A, 0)))
    tsub = A.period/A.nperiod
    tt, = quadgk(t -> tr(tpmeval(A,t)), 0., tsub; rtol)
    return tt/tsub
end
function LinearAlgebra.opnorm(A::FourierFunctionMatrix, p::Union{Real, Missing} = missing) 
    return convert(FourierFunctionMatrix,opnorm(convert(PeriodicFunctionMatrix,A),p))
end
function LinearAlgebra.norm(A::FourierFunctionMatrix, p::Real = 2; rtol = sqrt(eps())) 
    isconstant(A) && (return norm(tpmeval(A,0)))
    tsub = A.period/A.nperiod
    if p == 2
       nrm, = quadgk(t -> norm(tpmeval(A,t))^2, 0., tsub; rtol)
       return sqrt(nrm*A.nperiod)
    elseif isinf(p)
        return -optimize(t->-norm(tpmeval(A,t)),0.,tsub,Optim.Brent(),rel_tol = rtol).minimum
    elseif p == 1    
        nrm, = quadgk(t -> norm(tpmeval(A,t)), 0., tsub; rtol)
        return nrm*A.nperiod
    else
        throw(ArgumentError("only p-norms for p = 1, 2, or Inf are supported"))
    end
end
# function LinearAlgebra.norm(A::FourierFunctionMatrix, p::Real = 2; K = 128) 
#     isconstant(A) && (return norm(tpmeval(A, 0)))
#     nrm = zero(eltype(A))
#     Δ = A.period/A.nperiod/K
#     ts = Δ/2
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
function +(A::FourierFunctionMatrix, B::FourierFunctionMatrix)
    period = promote_period(A, B)
    #nperiod = numerator(rationalize(period/A.period))*A.nperiod
    domain(A.M) == domain(B.M) && (return FourierFunctionMatrix(A.M+B.M, period))
    FourierFunctionMatrix(Fun(t-> A.M(t),Fourier(0..period)),period)+FourierFunctionMatrix(Fun(t-> B.M(t),Fourier(0..period)),period)
end
+(A::FourierFunctionMatrix, C::AbstractMatrix) = +(A, FourierFunctionMatrix(C, A.period))
+(A::AbstractMatrix, C::FourierFunctionMatrix) = +(FourierFunctionMatrix(A, C.period), C)
+(A::FourierFunctionMatrix, C::PeriodicFunctionMatrix) = +(convert(PeriodicFunctionMatrix,A), C)
+(A::PeriodicFunctionMatrix, C::FourierFunctionMatrix) = +(A, convert(PeriodicFunctionMatrix,C))
-(A::FourierFunctionMatrix) = FourierFunctionMatrix(-A.M, A.period)
-(A::FourierFunctionMatrix, B::FourierFunctionMatrix) = +(A,-B)
-(A::FourierFunctionMatrix, C::AbstractMatrix) = +(A,-C)
-(A::AbstractMatrix, C::FourierFunctionMatrix) = +(A, -C)
function (+)(A::FourierFunctionMatrix, J::UniformScaling{<:Real}) 
    m, n = size(A)
    n == m || throw(DimensionMismatch("matrix is not square: dimensions are $((m,n))"))
    A+Matrix(J(n))
end
(+)(J::UniformScaling{<:Real}, A::FourierFunctionMatrix) = +(A,J)
(-)(A::FourierFunctionMatrix, J::UniformScaling{<:Real}) = +(A,-J)
(-)(J::UniformScaling{<:Real}, A::FourierFunctionMatrix) = +(-A,J)
-(A::FourierFunctionMatrix, C::PeriodicFunctionMatrix) = +(A,-C)
-(A::PeriodicFunctionMatrix, C::FourierFunctionMatrix) = +(A, -C)


function *(A::FourierFunctionMatrix, B::FourierFunctionMatrix)
    period = promote_period(A, B)
    domain(A.M) == domain(B.M)  && (return FourierFunctionMatrix(A.M*B.M, period))
    convert(FourierFunctionMatrix,convert(PeriodicFunctionMatrix,A) * convert(PeriodicFunctionMatrix,B))
end
*(A::FourierFunctionMatrix, C::AbstractMatrix) = *(A, FourierFunctionMatrix(C, A.period))
*(A::AbstractMatrix, C::FourierFunctionMatrix) = *(FourierFunctionMatrix(A, C.period), C)
*(A::FourierFunctionMatrix, C::Real) = FourierFunctionMatrix(C*A.M, A.period)
*(C::Real, A::FourierFunctionMatrix) = FourierFunctionMatrix(C*A.M, A.period)
/(A::FourierFunctionMatrix, C::Real) = *(A, 1/C)
*(J::UniformScaling{<:Real}, A::FourierFunctionMatrix) = J.λ*A
*(A::FourierFunctionMatrix, J::UniformScaling{<:Real}) = A*J.λ
*(A::FourierFunctionMatrix, C::PeriodicFunctionMatrix) = *(convert(PeriodicFunctionMatrix,A), C)
*(A::PeriodicFunctionMatrix, C::FourierFunctionMatrix) = *(A, convert(PeriodicFunctionMatrix,C))


function horzcat(A::FourierFunctionMatrix, B::FourierFunctionMatrix)
    A.period == B.period && A.nperiod == B.nperiod && (return FourierFunctionMatrix(Fun(t->[A.M(t) B.M(t)],Fourier(0..A.period)), A.period))
    convert(FourierFunctionMatrix,[convert(PeriodicFunctionMatrix,A) convert(PeriodicFunctionMatrix,B)])
end
hcat(A::FourierFunctionMatrix, B::FourierFunctionMatrix) = horzcat(A,B)
hcat(A::FourierFunctionMatrix, C::AbstractMatrix) = horzcat(A, FourierFunctionMatrix(C, A.period))
hcat(A::AbstractMatrix, C::FourierFunctionMatrix) = horzcat(FourierFunctionMatrix(A, C.period), C)
horzcat(A::FourierFunctionMatrix, C::AbstractMatrix) = horzcat(A, FourierFunctionMatrix(C, A.period))
horzcat(A::AbstractMatrix, C::FourierFunctionMatrix) = horzcat(FourierFunctionMatrix(A, C.period), C)


function vertcat(A::FourierFunctionMatrix, B::FourierFunctionMatrix)
    A.period == B.period && A.nperiod == B.nperiod && (return FourierFunctionMatrix(Fun(t->[A.M(t); B.M(t)],Fourier(0..A.period)), A.period))
    convert(FourierFunctionMatrix,[convert(PeriodicFunctionMatrix,A); convert(PeriodicFunctionMatrix,B)])
end
vcat(A::FourierFunctionMatrix, B::FourierFunctionMatrix) = vertcat(A,B)
vcat(A::FourierFunctionMatrix, C::AbstractMatrix) = vertcat(A, FourierFunctionMatrix(C, A.period))
vcat(A::AbstractMatrix, C::FourierFunctionMatrix) = vertcat(FourierFunctionMatrix(A, C.period), C)
vertcat(A::FourierFunctionMatrix, C::AbstractMatrix) = vertcat(A, FourierFunctionMatrix(C, A.period))
vertcat(A::AbstractMatrix, C::FourierFunctionMatrix) = vertcat(FourierFunctionMatrix(A, C.period), C)


function blockdiag(A::FourierFunctionMatrix, B::FourierFunctionMatrix)
    A.period == B.period && A.nperiod == B.nperiod && (return FourierFunctionMatrix(Fun(t->bldiag(A.M(t), B.M(t)),Fourier(0..A.period)), A.period))
    convert(FourierFunctionMatrix,blockdiag(convert(PeriodicFunctionMatrix,A), convert(PeriodicFunctionMatrix,B)))
end

Base.iszero(A::FourierFunctionMatrix) = all(iszero,coefficients(A.M))
LinearAlgebra.issymmetric(A::FourierFunctionMatrix) = iszero(A.M-transpose(A.M))
function ==(A::FourierFunctionMatrix, B::FourierFunctionMatrix)
    isconstant(A) && isconstant(B) && (return iszero(A.M(0)-B.M(0)))
    A.period*B.nperiod ≈ B.period*A.nperiod && domain(A.M) == domain(B.M) && iszero(A-B) #iszero(t->A.M(t)-B.M(t), A.period/A.nperiod) 
end
function Base.isapprox(A::FourierFunctionMatrix, B::FourierFunctionMatrix; rtol::Real = sqrt(eps(Float64)), atol::Real = 0)
    isconstant(A) && isconstant(B) && (return isapprox(A.M(0), B.M(0); rtol, atol))
    ts = rand()*A.period
    A.period*B.nperiod ≈ B.period*A.nperiod && domain(A.M) == domain(B.M) && isapprox(A.M(ts), B.M(ts); rtol, atol) 
end
function Base.isapprox(A::FourierFunctionMatrix, J::UniformScaling{<:Real}; kwargs...)
    isconstant(A) && (return isapprox(A.M(0), J; kwargs...))
    ts = rand()*A.period
    isapprox(tpmeval(A,ts), J; kwargs...) 
end
Base.isapprox(J::UniformScaling{<:Real}, A::FourierFunctionMatrix; kwargs...) = isapprox(A, J; kwargs...)


