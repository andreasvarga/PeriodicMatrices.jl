# eye(n) = Matrix{Bool}(I, n, n)
# eye(m,n) = Matrix{Bool}(I, m, n)
# eye(::Type{T}, n) where {T} = Matrix{T}(I, n, n)
eye(::Type{T}, m, n) where {T} = Matrix{T}(I, m, n)
function pmzeros(::Type{T},m::Vector{Int},n::Vector{Int}) where {T}
    lm = length(m)
    ln = length(n)
    return [zeros(T,m[mod(i-1,lm)+1], n[mod(i-1,ln)+1]) for i in 1:lcm(lm,ln)]
end
pmzeros(m::Vector{Int},n::Vector{Int}) = pmzeros(Float64,m,n)


"""
     ev = pseig(A::Array{T,3}; rev = true, fast = false) 

Compute the eigenvalues of a product of `p` square matrices 
`A(p)...*A(2)*A(1)`, if `rev = true` (default) (also called characteristic multipliers) or 
of `A(1)*A(2)...A(p)` if `rev = false`, without evaluating the product. 
The matrices `A(1)`, `...`, `A(p)` are contained in the `n×n×p` array `A` 
such that the `i`-th matrix `A(i)` is contained in `A[:,:,i]`.
If `fast = false` (default) then the eigenvalues are computed using an approach
based on the periodic Schur decomposition [1], while if `fast = true` 
the structure exploiting reduction [2] of an appropriate lifted pencil is employed.
This later option may occasionally lead to inaccurate results for large number of matrices. 

_References_

[1] A. Bojanczyk, G. Golub, and P. Van Dooren, 
    The periodic Schur decomposition. Algorithms and applications, Proc. SPIE 1996.

[2] A. Varga & P. Van Dooren. Computing the zeros of periodic descriptor systems.
    Systems and Control Letters, 50:371-381, 2003.

"""
function pseig(A::Array{T,3}; rev::Bool = true, fast::Bool = false) where T
   n = size(A,1)
   n == size(A,2) || error("A must have equal first and second dimensions") 
   p = size(A,3)
   if fast 
      if rev 
         ev = eigvals(psreduc_reg(A)...)
      else
         imap = p:-1:1                     
         ev = eigvals(psreduc_reg(view(A,:,:,imap))...)
      end
      isreal(ev) && (ev = real(ev))
      sorteigvals!(ev)
      return sort!(ev,by=abs,rev=true)
   else
      T1 = promote_type(Float64,T)
      ev = pschur(T1.(A); rev, withZ = false)[3]
      isreal(ev) && (ev = real(ev))
      return ev
   end
end
"""
     ev = pseig(A::Vector{Matrix}[, k = 1]; rev = true, fast = false) 

Compute the eigenvalues of a square cyclic product of `p` matrices 
`A(k-1)...*A(2)*A(1)*A(p)...*A(k)`, if `rev = true` (default) or 
`A(k)*A(k+1)*...A(p)*A(1)...A(k-1)` if `rev = false`, without evaluating the product. 
The argument `k` specifies the starting index (default: `k = 1`). 
The matrices `A(1)`, `...`, `A(p)` are contained in the `p`-vector of matrices `A` 
such that the `i`-th matrix  `A(i)`, of dimensions `m(i)×n(i)`, is contained in `A[i]`.
If `fast = false` (default) then the eigenvalues are computed using an approach
based on the periodic Schur decomposition [1], while if `fast = true` 
the structure exploiting reduction [2] of an appropriate lifted pencil is employed. 
This later option may occasionally lead to inaccurate results for large number of matrices. 

_Note:_ The first `nmin` components of `ev` contains the `core eigenvalues` of the appropriate matrix product,
where `nmin` is the minimum row dimensions of matrices `A[i]`, for `i = 1, ..., p`, 
while the last `ncur-nmin` components of `ev` are zero, 
where `ncur` is the column dimension of `A[k]` if `rev = true` or 
the row dimension of `A[k]` if `rev = false`. 

_References_

[1] A. Bojanczyk, G. Golub, and P. Van Dooren, 
    The periodic Schur decomposition. Algorithms and applications, Proc. SPIE 1996.

[2] A. Varga & P. Van Dooren. Computing the zeros of periodic descriptor systems.
    Systems and Control Letters, 50:371-381, 2003.

"""
function pseig(A::Vector{Matrix{T}}, k::Int = 1; rev::Bool = true, fast::Bool = false) where T
   p = length(A)
   istart = mod(k-1,p)+1
   nev = rev ? size(A[istart],2) : size(A[istart],1)
   # check dimensions
   m, n = size.(A,1), size.(A,2)
   if rev
      all(m .== view(n,mod.(1:p,p).+1)) || error("incompatible dimensions")
   else
      all(n .== view(m,mod.(1:p,p).+1)) || error("incompatible dimensions")
   end
   if fast 
      ncore = minimum(size.(A,1))
      if istart == 1 && rev 
         ev = eigvals(psreduc_reg(A)...)
      else
         imap = rev ? (mod.(istart:istart+p-1,p).+1) :
                      (mod.(p-istart+1:-1:-istart+2,p).+1)                     
         ev = eigvals(psreduc_reg(view(A,imap))...)
      end
      isreal(ev) && (ev = real(ev))
      T <: Complex || sorteigvals!(ev)
      ind = sortperm(ev; by = abs, rev = true) # select the core eigenvalues
      return [ev[ind[1:ncore]]; zeros(eltype(ev),nev-ncore)]  # pad with the necessary zeros
   else
      if istart == 1 
         ev = pschur(A; rev, withZ = false)[3]
      else
         # avoid repeated reindexing
         imap = mod.(istart-1:istart+p-2,p).+1
         rev && reverse!(imap)        
         ev = pschur(view(A,imap); rev = false, withZ = false)[3]
      end
      isreal(ev) && (ev = real(ev))
      return ev[1:nev]
   end
end

function sorteigvals!(ev)
   # an approximately complex conjugated set is assumed 
   isreal(ev) && (return ev)
   tc = ev[imag.(ev) .> 0]
   ev[:] = [ev[imag.(ev) .== 0]; sort([tc; conj.(tc)],by = real)]
   return ev
end


# conversions
"""
     ts2fm(A::Vector{<:AbstractMatrix}, period; method = "linear") -> At::Function

Compute the function matrix corresponding to an interpolated matrix time series. 
For the given matrix time series `A`, a function matrix `A(t)` is defined as the 
mapping `A(t) = t -> etpf(t)`, where `etpf(t)` is an interpolation object,  
as provided in the [`Interpolations.jl`](https://github.com/JuliaMath/Interpolations.jl)  package. 
The keyword parameter `method` specifies the interpolation method to be used as follows:

`method = "constant"` - use periodic B-splines of degree 0 (constant interpolation);

`method = "linear"` - use periodic B-splines of degree 1 (linear interpolation) (default);

`method = "quadratic"` - use periodic B-splines of degree 2 (quadratic interpolation); 

`method = "cubic"` - use periodic B-splines of degree 3 (cubic interpolation). 
"""
function ts2fm(A::Vector{<:AbstractMatrix}, T; method = "linear")
   N = length(A)
   N == 0 && error("empty time array not supported")
   N == 1 && (return t -> A[1])
   #ts = (0:N-1)*(A.period/N)
   d = T/(N-1)
   ts = (0:N-1)*d
   n1, n2 = size(A[1])
   intparray = Array{Any,2}(undef, n1, n2)
   if method == "constant"     
      # use simple function evaluation
      [intparray[i,j] = scale(interpolate(getindex.(A,i,j), BSpline(Constant())), ts) for i in 1:n1, j in 1:n2]
   elseif method == "linear" || N == 2
      [intparray[i,j] = scale(interpolate(getindex.(A,i,j), BSpline(Linear())), ts) for i in 1:n1, j in 1:n2]
   elseif method == "quadratic" || N == 3     
      [intparray[i,j] = scale(interpolate(getindex.(A,i,j), BSpline(Quadratic())), ts) for i in 1:n1, j in 1:n2]
   elseif method == "cubic"     
      [intparray[i,j] = scale(interpolate(getindex.(A,i,j), BSpline(Cubic())), ts) for i in 1:n1, j in 1:n2]
   else
      throw(ArgumentError("no such option method = $method"))
   end
   return t -> [intparray[i,j](t) for i in 1:n1, j in 1:n2 ]
end

function psreduc_fast(S::Vector{Matrix{T1}}, T::Vector{Matrix{T1}}; atol::Real = 0) where T1
    # PSREDUC_FAST  Finds a matrix pair having the same finite and infinite 
    #               eigenvalues as a given periodic pair.
    #
    # [A,E] = PSREDUC_FAST(S,T,tol) computes a matrix pair (A,E) having 
    # the same finite and infinite eigenvalues as the periodic pair (S,T).
    #
    
    #   A. Varga 30-3-2004. 
    #   Revised: .
    #
    #   Reference:
    #   A. Varga & P. Van Dooren
    #      Computing the zeros of periodic descriptor systems.
    #      Systems and Control Letters, vol. 50, 2003.
    
    K = length(S) 
    if K == 1
       return S[1], T[1] 
    else
       m = size.(S,1)
       n = size.(S,2)
       if sum(m) > sum(n) 
           # build the dual pair
        #    [S,T]=celltranspose(S,T); 
        #    [S{1:K}] = deal(S{K:-1:1});
        #    [T{1:K-1}] = deal(T{K-1:-1:1});
        #    m = size.(S,1)
        #    n = size.(S,2)
       end 
    
       si = S[1];  ti = -T[1];
       tolr = atol
       for i = 1:K-2
           F = qr([ ti; S[i+1] ], ColumnNorm()) 
           nr = minimum(size(F.R))
           # compute rank of r 
           ss = abs.(diag(F.R[1:nr,1:nr]))
           atol == 0 && ( tolr = (m[i]+m[i+1]) * maximum(ss) * eps())
           rankr = count(ss .> tolr)
    
           si = F.Q'*[si; zeros(T1,m[i+1],n[1])]; si=si[rankr+1:end,:]
           ti = F.Q'*[ zeros(T1,size(ti,1),n[i+2]); -T[i+1]]; ti = ti[rankr+1:end,:]
       end
       a = [ zeros(T1,m[K],n[1]) S[K]; si ti] 
       e = [ T[K] zeros(T1,m[K],n[K]); zeros(size(si)...) zeros(size(ti)...)] 
       return a, e
    end
end
"""
    psreduc_reg(A) -> (M, N)

Determine for a `n×n×p` array `A`, the matrix pair `(M, N)` 
with `N` invertible and `M-λN` regular, such that 
the eigenvalues of `M-λN` are the same as those of the matrix product
`A(p)*A(p-1)*...*A(1)`, where `A(i)` is contained in `A[:,:,i]`. 
The structure exploiting fast reduction method of [1] is employed to determine `M` and `N`.

[1] A. Varga & P. Van Dooren. Computing the zeros of periodic descriptor systems.
    Systems and Control Letters, 50:371-381, 2003.
"""
function psreduc_reg(A::AbstractArray{T,3}) where T
     
    K = size(A,3) 
    n = size(A,1)
    n == size(A,2) || error("A must have equal first and second dimensions")
    if K == 1
       return A[:,:,1], Matrix{T}(I, n, n)
    else   
       Z = zeros(T,n,n)
       ZI = [Z; -I]
       si = A[:,:,1];  ti = -I
       for i = 1:K-1
           F = qr([ ti; A[:,:,i+1] ])     
           si = F.Q'*[si; Z];  si = si[n+1:end,:]
           ti = F.Q'*ZI; ti = ti[n+1:end,:]
       end
       return si, -ti
    end
end
"""
    psreduc_reg(A,E) -> (M, N)

Determine for a pair of `n×n×p` arrays `(A,E)`, the matrix pair `(M, N)` 
with `M-λN` regular, such that 
the eigenvalues of `M-λN` are the same as those of the quotient matrix product
`inv(E(p))*(A(p)*inv(E(p-1))*A(p-1)*...*inv(E(1))*A(1)`, where `A(i)` is contained in `A[:,:,i]` and `E(i)` is contained in `E[:,:,i]`. 
The structure exploiting fast reduction method of [1] is employed to determine `M` and `N`.

[1] A. Varga & P. Van Dooren. Computing the zeros of periodic descriptor systems.
    Systems and Control Letters, 50:371-381, 2003.
"""
function psreduc_reg(A::AbstractArray{T,3}, E::AbstractArray{T,3}) where T
     
    K = size(A,3) 
    n = size(A,1)
    n == size(A,2) || throw(ArgumentError("A must have equal first and second dimensions"))
    (n, n, K) == size(E) || throw(ArgumentError("A and E must have the same dimensions"))

    if K == 1
       return A[:,:,1], E[:,:,1] 
    else   
       Z = zeros(T,n,n)
       si = A[:,:,1];  ti = -E[:,:,1]
       for i = 1:K-1
           F = qr([ ti ; A[:,:,i+1] ])     
           si = F.Q'*[si; Z];  si = si[n+1:end,:]
           ti = F.Q'*[Z; -E[:,:,i+1]]; ti = ti[n+1:end,:]
       end
       return si, -ti
    end
end

function psreduc_reg(A::AbstractVector{Matrix{T}}, E::AbstractVector{Matrix{T}}) where {T}
     
   K = length(A) 
   K == length(E) || throw(ArgumentError("A and E must have the same lengths"))
   nd = size.(A,1); n = size.(A,2)
   # n == nd[mod.(-1:K-2,K).+1] || 
   #    error("the number of columns of A[i] must be equal to the number of rows of A[i-1]")
   # nde = size.(E,1); ne = size.(E,2)
   # all(nde .== ne) || error("all E[i] must be square")
   # all(ne .== nd) || error("the number of rows of A[i] must be equal to the order of E[i]")

   if K == 1
      return A[1], E[1] 
   else   
      si = A[1];  ti = -E[1]
      for i = 1:K-1
          F = qr([ ti ; A[i+1] ])   
          ni1 = n[i+1] 
          si = F.Q'*[si; zeros(T,nd[i+1],n[1])];  si = si[ni1+1:end,:] 
          ip2 = i+2; ip2 > K && (ip2 = 1)
          ti = F.Q'*[zeros(T,size(ti,1),n[ip2]); -E[i+1]]; ti = ti[ni1+1:end,:] 
      end
      return si, -ti
   end
end


"""
    psreduc_reg(A) -> (M, N)

Determine for a `p`-dimensional vector of rectangular matrices `A`, 
the matrix pair `(M, N)` with `N` invertible and `M-λN` regular, such that 
the eigenvalues of `M-λN` are the same as those of the square 
matrix product `A(p)*A(p-1)*...*A(1)`, where `A(i)` is contained in `A[i]`. 
The structure exploiting fast reduction method of [1] is employed to determine `M` and `N`.

[1] A. Varga & P. Van Dooren. Computing the zeros of periodic descriptor systems.
    Systems and Control Letters, 50:371-381, 2003.
"""
function psreduc_reg(A::AbstractVector{Matrix{T}}) where T
     
   K = length(A) 
   n = size.(A,2) 
   n == size.(A,1)[mod.(-1:K-2,K).+1] || 
      error("the number of columns of A[i] must be equal to the number of rows of A[i-1]")
   k = mod(1,K)+1
   si = A[1];  ti = Matrix{T}(I, n[k], n[k])
   for i = 2:K
       F = qr([ -ti; A[i] ])  
       mi1 = n[mod(i,K)+1]
       si = (F.Q'*[si; zeros(T,mi1,n[1])])[n[i]+1:end,:] 
       ti = (F.Q'*[ zeros(T,n[i],mi1); I])[n[i]+1:end,:] 
   end
   return si, ti
end

