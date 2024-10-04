# Operations with periodic symbolic matrices

function pmderiv(A::PeriodicSymbolicMatrix) 
    @variables t   
    return PeriodicSymbolicMatrix{:c,Num}(Symbolics.derivative(A.F,t), A.period, nperiod = A.nperiod)
end
#LinearAlgebra.inv(A::PeriodicSymbolicMatrix) = PeriodicSymbolicMatrix(inv(A.F), A.period; nperiod = A.nperiod)
function LinearAlgebra.inv(A::PeriodicSymbolicMatrix)
    PeriodicSymbolicMatrix(inv(A.F), A.period; nperiod = A.nperiod)
end
function LinearAlgebra.transpose(A::PeriodicSymbolicMatrix)  
    return PeriodicSymbolicMatrix{:c,Num}(copy(transpose(A.F)), A.period, nperiod = A.nperiod)
end
function LinearAlgebra.adjoint(A::PeriodicSymbolicMatrix)  
    return PeriodicSymbolicMatrix{:c,Num}(copy(adjoint(A.F)), A.period, nperiod = A.nperiod)
end
function Symbolics.simplify(A::PeriodicSymbolicMatrix)
    return PeriodicSymbolicMatrix{:c,Num}(Symbolics.simplify.(A.F), A.period; nperiod = A.nperiod)
end
LinearAlgebra.tr(A::PeriodicSymbolicMatrix) = PeriodicSymbolicMatrix([tr(A.F)], A.period; nperiod = A.nperiod)
# function trace(A::PeriodicSymbolicMatrix; K = 128) 
#     @variables t 
#     trs = tr(A.F)  
#     isconstant(A) && (return Symbolics.unwrap.(substitute.(trs, (Dict(t => 0),)))[1]*Δ)
#     ts = zero(eltype(A))
#     Δ = A.period/A.nperiod/K
#     tt = zero(eltype(Δ))
#     for i = 1:K
#         tt += Symbolics.unwrap.(substitute.(trs, (Dict(t => ts),)))[1]*Δ
#         ts += Δ
#     end 
#     return tt*A.nperiod/A.period
# end
function trace(A::PeriodicSymbolicMatrix; rtol=sqrt(eps())) 
    isconstant(A) && (return tr(tpmeval(A, 0)))
    tsub = A.period/A.nperiod
    tt, = quadgk(t -> tr(tpmeval(A,t)), 0., tsub; rtol)
    return tt/tsub
end
function LinearAlgebra.opnorm(A::PeriodicSymbolicMatrix, p::Real = 2) 
    if p == 2
       return PeriodicSymbolicMatrix([norm(A.F,2)],A.period; nperiod = A.nperiod) 
    elseif p == 1 || isinf(p)
       return PeriodicSymbolicMatrix([opnorm(A.F,p)],A.period; nperiod = A.nperiod) 
    else
       throw(ArgumentError("only p-opnorms for p = 1, 2, or Inf are supported"))
    end
end

function LinearAlgebra.norm(A::PeriodicSymbolicMatrix, p::Real = 2; rtol = sqrt(eps())) 
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
# function LinearAlgebra.norm(A::PeriodicSymbolicMatrix, p::Real = 2; K = 128) 
#     @variables t   
#     isconstant(A) && (return norm(Symbolics.unwrap.(substitute.(A.F, (Dict(t => 0),)))))
#     nrm = zero(eltype(A))
#     Δ = A.period/A.nperiod/K
#     ts = zero(eltype(Δ))
#     if p == 2
#        for i = 1:K
#            nrm += norm(Symbolics.unwrap.(substitute.(A.F, (Dict(t => ts),))))^2*Δ
#            ts += Δ
#        end 
#        return sqrt(nrm*A.nperiod)
#     elseif isinf(p)
#         for i = 1:K
#             nrm = max(nrm,norm(Symbolics.unwrap.(substitute.(A.F, (Dict(t => ts),)))))
#             ts += Δ
#         end 
#         return nrm
#     elseif p == 1    
#         for i = 1:K
#             nrm += norm(Symbolics.unwrap.(substitute.(A.F, (Dict(t => ts),))))*Δ
#             ts += Δ
#         end 
#         return nrm*A.nperiod
#     else
#         throw(ArgumentError("only p-norms for p = 1, 2, or Inf are supported"))
#     end
# end 
function +(A::PeriodicSymbolicMatrix, B::PeriodicSymbolicMatrix)
    period, nperiod = promote_period2(A, B)
    # nperiod = numerator(rationalize(period/A.period))*A.nperiod
    #nperiod = gcd(A.nperiod,B.nperiod)
    #return PeriodicSymbolicMatrix{:c,Num}(Symbolics.simplify.(A.F + B.F), period; nperiod)
    return PeriodicSymbolicMatrix{:c,Num}(A.F + B.F, period; nperiod)
end
+(A::PeriodicSymbolicMatrix, C::AbstractMatrix) = +(A, PeriodicSymbolicMatrix(C, A.period))
+(A::AbstractMatrix, C::PeriodicSymbolicMatrix) = +(PeriodicSymbolicMatrix(A, C.period), C)
+(A::PeriodicSymbolicMatrix, C::PeriodicFunctionMatrix) = +(convert(PeriodicFunctionMatrix,A), C)
+(A::PeriodicFunctionMatrix, C::PeriodicSymbolicMatrix) = +(A, convert(PeriodicFunctionMatrix,C))
-(A::PeriodicSymbolicMatrix) = PeriodicSymbolicMatrix(-A.F, A.period; nperiod = A.nperiod)
-(A::PeriodicSymbolicMatrix, B::PeriodicSymbolicMatrix) = +(A,-B)
-(A::PeriodicSymbolicMatrix, C::AbstractMatrix) = +(A,-C)
-(A::AbstractMatrix, C::PeriodicSymbolicMatrix) = +(A, -C)
-(A::PeriodicSymbolicMatrix, C::PeriodicFunctionMatrix) = +(A,-C)
-(A::PeriodicFunctionMatrix, C::PeriodicSymbolicMatrix) = +(A, -C)

function (+)(A::PeriodicSymbolicMatrix, J::UniformScaling{<:Real}) 
    m, n = size(A)
    n == m || throw(DimensionMismatch("matrix is not square: dimensions are $((m,n))"))
    PeriodicSymbolicMatrix(A.F+Matrix(J(n)), A.period; nperiod = A.nperiod)
end
(+)(J::UniformScaling{<:Real}, A::PeriodicSymbolicMatrix) = +(A,J)
(-)(A::PeriodicSymbolicMatrix, J::UniformScaling{<:Real}) = +(A,-J)
(-)(J::UniformScaling{<:Real}, A::PeriodicSymbolicMatrix) = +(-A,J)

function *(A::PeriodicSymbolicMatrix, B::PeriodicSymbolicMatrix)
    period, nperiod = promote_period2(A, B)
    # nperiod = numerator(rationalize(period/A.period))*A.nperiod
    return PeriodicSymbolicMatrix{:c,Num}(A.F * B.F, period; nperiod)
 end
*(A::PeriodicSymbolicMatrix, C::AbstractMatrix) = *(A, PeriodicSymbolicMatrix(C, A.period))
*(A::AbstractMatrix, C::PeriodicSymbolicMatrix) = *(PeriodicSymbolicMatrix(A, C.period), C)
*(A::PeriodicSymbolicMatrix, C::PeriodicFunctionMatrix) = *(convert(PeriodicFunctionMatrix,A), C)
*(A::PeriodicFunctionMatrix, C::PeriodicSymbolicMatrix) = *(A, convert(PeriodicFunctionMatrix,C))
*(A::PeriodicSymbolicMatrix, C::Real) = PeriodicSymbolicMatrix(C*A.F, A.period; nperiod = A.nperiod)
*(C::Real, A::PeriodicSymbolicMatrix) = PeriodicSymbolicMatrix(C*A.F, A.period; nperiod = A.nperiod)
/(A::PeriodicSymbolicMatrix, C::Real) = *(A, 1/C)
*(J::UniformScaling{<:Real}, A::PeriodicSymbolicMatrix) = J.λ*A 
*(A::PeriodicSymbolicMatrix, J::UniformScaling{<:Real}) = A*J.λ

function horzcat(A::PeriodicSymbolicMatrix, B::PeriodicSymbolicMatrix)
    period, nperiod = promote_period2(A, B)
    # nperiod = numerator(rationalize(period/A.period))*A.nperiod
    return PeriodicSymbolicMatrix{:c,Num}([A.F B.F], period; nperiod)
end
hcat(A::PeriodicSymbolicMatrix, B::PeriodicSymbolicMatrix) = horzcat(A,B)
hcat(A::PeriodicSymbolicMatrix, C::AbstractMatrix) = horzcat(A, PeriodicSymbolicMatrix(C, A.period))
hcat(A::AbstractMatrix, C::PeriodicSymbolicMatrix) = horzcat(PeriodicSymbolicMatrix(A, C.period), C)
horzcat(A::PeriodicSymbolicMatrix, C::AbstractMatrix) = horzcat(A, PeriodicSymbolicMatrix(C, A.period))
horzcat(A::AbstractMatrix, C::PeriodicSymbolicMatrix) = horzcat(PeriodicSymbolicMatrix(A, C.period), C)


function vertcat(A::PeriodicSymbolicMatrix, B::PeriodicSymbolicMatrix)
    period, nperiod = promote_period2(A, B)
    # nperiod = numerator(rationalize(period/A.period))*A.nperiod
    return PeriodicSymbolicMatrix{:c,Num}([A.F; B.F], period; nperiod)
end
vcat(A::PeriodicSymbolicMatrix, B::PeriodicSymbolicMatrix) = vertcat(A,B)
vcat(A::PeriodicSymbolicMatrix, C::AbstractMatrix) = vertcat(A, PeriodicSymbolicMatrix(C, A.period))
vcat(A::AbstractMatrix, C::PeriodicSymbolicMatrix) = vertcat(PeriodicSymbolicMatrix(A, C.period), C)
vertcat(A::PeriodicSymbolicMatrix, C::AbstractMatrix) = vertcat(A, PeriodicSymbolicMatrix(C, A.period))
vertcat(A::AbstractMatrix, C::PeriodicSymbolicMatrix) = vertcat(PeriodicSymbolicMatrix(A, C.period), C)


function blockdiag(A::PeriodicSymbolicMatrix, B::PeriodicSymbolicMatrix)
    period, nperiod = promote_period2(A, B)
    # nperiod = numerator(rationalize(period/A.period))*A.nperiod
    return PeriodicSymbolicMatrix{:c,Num}(bldiag(A.F, B.F), period; nperiod)
end



Base.iszero(A::PeriodicSymbolicMatrix) = iszero(A.F)
LinearAlgebra.issymmetric(A::PeriodicSymbolicMatrix) = iszero(A.F-transpose(A.F))
function ==(A::PeriodicSymbolicMatrix, B::PeriodicSymbolicMatrix)
    isconstant(A) && isconstant(B) && (return iszero(A.F-B.F))
    A.period*B.nperiod ≈ B.period*A.nperiod && iszero(A.F-B.F) 
end
==(A::PeriodicSymbolicMatrix, J::UniformScaling{<:Real}) = iszero(A.F-J) 
==(J::UniformScaling{<:Real},A::PeriodicSymbolicMatrix) = iszero(A.F-J) 
function Base.isapprox(A::PeriodicSymbolicMatrix, B::PeriodicSymbolicMatrix; rtol::Real = sqrt(eps(Float64)), atol::Real = 0)
    A.period*B.nperiod ≈ B.period*A.nperiod || (return false)
    @variables t   
    ts = rand()*A.period
    isconstant(A) && isconstant(B) && (return isapprox(Symbolics.unwrap.(substitute.(A.F, (Dict(t => ts),))), Symbolics.unwrap.(substitute.(B.F, (Dict(t => ts),))); rtol, atol))
    isapprox(Symbolics.unwrap.(substitute.(A.F, (Dict(t => ts),))), Symbolics.unwrap.(substitute.(B.F, (Dict(t => ts),))); rtol, atol) 
end
function Base.isapprox(A::PeriodicSymbolicMatrix, J::UniformScaling{<:Real}; kwargs...)
    @variables t   
    ts = rand()*A.period
    isapprox(Symbolics.unwrap.(substitute.(A.F, (Dict(t => ts),))), J; kwargs...) 
end
Base.isapprox(J::UniformScaling{<:Real}, A::PeriodicSymbolicMatrix; kwargs...) = isapprox(A, J; kwargs...)
