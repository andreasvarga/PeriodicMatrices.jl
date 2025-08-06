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
     ev = peigvals(A::Array{T,3}; rev = true, fast = false) 

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
function peigvals(A::Array{T,3}; rev::Bool = true, fast::Bool = false) where T
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
     ev = peigvals(A::Vector{Matrix}[, k = 1]; rev = true, fast = false) 

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
function peigvals(A::Vector{Matrix{T}}, k::Int = 1; rev::Bool = true, fast::Bool = false, PSD_SLICOT::Bool = true) where T
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
      if PSD_SLICOT
         if istart == 1 
            ev = pschur(A; rev, withZ = false)[3]
         else
            # avoid repeated reindexing
            imap = mod.(istart-1:istart+p-2,p).+1
            rev && reverse!(imap)        
            ev = pschur(view(A,imap); rev = false, withZ = false)[3]
         end
      else
         if istart == 1 
            ev = PeriodicSchurDecompositions.pschur(A, rev ? (:L) : (:R); wantZ = false, wantT = false).values
         else
            # avoid repeated reindexing
            imap = mod.(istart-1:istart+p-2,p).+1
            rev && reverse!(imap)        
            ev = PeriodicSchurDecompositions.pschur(view(A,imap), :L ; wantZ = false, wantT = false).values
         end
      end
      isreal(ev) && (ev = real(ev))
      return ev[1:nev]
   end
end

"""
    peigvecs(A::Vector{Matrix}; rev = true, select_fun = f(x), allvecs::Bool) => (V::Vector{Matrix}, ev::Vector{Complex})

Compute eigenvalues `ev` and the corresponding right eigenvectors `V` of a cyclic product of `p` `n×n` real matrices 
`A(p)*...*A(2)*A(1) =: Ψ`, if `rev = true` (default) or 
`A(1)*A(2)*...*A(p) =: Ψ` if `rev = false`, without evaluating the product. 
The vectors are returned as columns of the elements of the vector `V` of matrices.    
The resulting vectors satisfy `A[k]*V[k] = V[k+1]*Γ[k]` if `rev = true`, or 
`A[k]*V[k+1] = V[k]*Γ[k]` if `rev = false`, where `Γ[k]` are diagonal matrices and satisfy `Γ[1]*Γ[2]...Γ[p] = Diagonal(ev)`
If `allvecs = false`, then `V` is a vector with a single component which satisfies
`Ψ*V[1] = V[1]*Diagonal(ev)`. 

A selection of eigenvectors can be computed corresponding to eigenvalues which satisfy `f(x) = true`, where `x` is an eigenvalue and `f(x)` is
a function of scalar parameter `x` (default: `f(x) = true`, thus all eigenvectors are computed).    

The eigenvalues are computed using an approach based on the periodic Schur decomposition [1], for which purpose SLICOT based
wrappers are used if `PSD_SLICOT = true` (default), or the generic software from the [`PeriodicSchurDecomposition.jl`](https://github.com/RalphAS/PeriodicSchurDecompositions.jl) package is employed, 
if `PSD_SLICOT = false`.
The eigenvectors are computed using the [`eigvecs`](https://github.com/RalphAS/PeriodicSchurDecompositions.jl/blob/676582eb0c0f1e9d9260f1852b841a3e090aff65/src/vectors.jl#L2) 
function available in the [`PeriodicSchurDecomposition.jl`](https://github.com/RalphAS/PeriodicSchurDecompositions.jl) package.    

_References_

[1] A. Bojanczyk, G. Golub, and P. Van Dooren, 
    The periodic Schur decomposition. Algorithms and applications, Proc. SPIE 1996.
"""
function peigvecs(A::Vector{Matrix{T}}, k::Int = 1; rev::Bool = true, PSD_SLICOT::Bool = true, select_fun = x->true, allvecs = true) where {T <: Real}
   mp, np = size.(A,1), size.(A,2) 
   n = maximum(np); m = maximum(mp)
   (n == minimum(np) && m == minimum(mp)) || throw(ArgumentError("only constant dimensions are supported"))
   n == m || throw(ArgumentError("all component matrices must be square"))
   if PSD_SLICOT
      S, Z, ev, ischur,  = pschur(A; rev, withZ = true)
      PSF = PeriodicSchur(S[ischur],S[[1:ischur-1;ischur+1:length(S)]], Z, ev, rev ? 'L' : 'R', ischur)
   else
      PSF = PeriodicSchurDecompositions.pschur(A, rev ? (:L) : (:R))
   end
   select = select_fun.(PSF.values)
   vecs = _eigvecs(PSF, select; shifted = allvecs)
   return vecs, PSF.values[select]
end




function sorteigvals!(ev)
   # an approximately complex conjugated set is assumed 
   isreal(ev) && (return ev)
   tc = ev[imag.(ev) .> 0]
   ev[:] = [ev[imag.(ev) .== 0]; sort([tc; conj.(tc)],by = real)]
   return ev
end

"""
    _eigvecs(ps::PeriodicSchur, select::Vector{Bool}; shifted::Bool) => V::Vector{Matrix}

Compute selected right eigenvectors of a product of matrices in a periodic Schur form.
This is a modified version of the `eigvecs` function from the 
[`PeriodicSchurDecompositions.jl`](https://github.com/RalphAS/PeriodicSchurDecompositions.jl) package.    

For the `n×n` matrices `Aⱼ, j ∈ 1:p`, whose periodic Schur decomposition
is in `ps`, this function computes vectors `v` such that
`Aₚ*...*A₂*A₁*v = λₖ v`  (or `A₁*A₂*...*Aₚ*v = λₖ v`) for left (right)
oriented `ps`. The returned vectors `V` correspond to the 
selected eigenvalues from `ps.values` via the `true` entries of `select`.

If keyword `shifted` is true (the default), eigenvectors for circularly shifted
permutations of the `A` matrices are also returned.
The vectors are returned as columns of the elements of the vector `V`
of matrices. (*Note:* A vector of length one is returned if `shifted` is false.)
The vectors are satisfy `Aⱼvⱼ = μⱼvⱼ₊₁`, for a left oriened product, or `Aⱼvⱼ₊₁ = μⱼvⱼ`
for a right oriented product, where `μ₁μ₂...μₚ = λₖ`.

If the element type of `ps` is real, `select` may be updated to include
conjugate pairs. `ps` itself is not modified.

The algorithm may fail if some selected eigenvalues are associated
with an invariant subspace that cannot be untangled.
"""
function _eigvecs(ps0::PeriodicSchur{T}, select::AbstractVector{Bool};
                 shifted=true, verbosity=1) where {T}
    RT = real(T)
    CT = complex(T)
    ps = deepcopy(ps0)
    if (length(ps.Z) == 0) || (size(ps.Z[1], 1) == 0)
        throw(ArgumentError("eigvecs requires Schur vectors in the PSD"))
    end
    n,m = size(ps.Z[1])
    if length(select) != m
        throw(ArgumentError("length of `select` must correspond to rank of Schur (sub-)space"))
    end
    p = ps.period
    left = ps.orientation == 'L'
    if !all(select)
        if T <: Real
            inpair = false
            for j in 1:m
                if inpair
                    if select[j - 1]
                        if verbosity > 0 && !select[j]
                            @info "adding $j to select for conjugacy"
                        end
                        select[j] = true
                    elseif select[j]
                        if verbosity > 0 && !select[j - 1]
                            @info "adding $(j-1) to select for conjugacy"
                        end
                        select[j - 1] = true
                    end
                    inpair = false
                    continue
                end
                inpair = !isreal(ps.values[j])
            end
        end
        ordschur!(ps, select)
    end
    nvec = count(select)
    sel = falses(m)
    sel[1:nvec] .= true
    nmat = shifted ? p : 1
    Vs = [Matrix{CT}(undef, n, nvec) for _ in 1:nmat]
    iλ = 1
    while iλ <= nvec
        if T <: Real && !isreal(ps.values[1])
            # I have a hammer so this must be a nail.
            # set up and solve the 2x2 cyclic problem
            μ = (ps.values[1] + 0im) ^ (1/RT(p))
            Zd = [-μ * Matrix{CT}(I, 2, 2) for _ in 1:p]
            Zl = [Matrix{CT}(undef, 2, 2) for _ in 1:p]
            il = 0
            for l in 1:p
                lx = left ? l : (p + 1 - l)
                if l == ps.schurindex
                    Tl = ps.T1[1:2,1:2]
                else
                    il += 1
                    Tl = ps.T[il][1:2,1:2]
                end
                Zl[lx] .= Tl
            end
            nsolve = 2
            rowx = 1
            colx = 1:nsolve
            y = zeros(CT, nsolve * p)
            # replace a row
            y[rowx] = 1
            Zd[1][rowx, :] .= 0
            Zl[p][rowx, :] .= 0
            Zd[1][rowx, colx] .= 1
            R, Zu, Zr, _ = PeriodicSchurDecompositions._babd_qr!(Zd, Zl, y)
            x = PeriodicSchurDecompositions._babd_solve!(R, Zu, Zr, y)
            t = 1 / norm(view(x, 1:nsolve))
            for l in 1:nmat
                if left
                    i0 = (l - 1) * nsolve
                else
                    i0 = l == 1 ? 0 : (p + 1 - l) * nsolve
                end
                vl = Vs[l]
                mul!(view(vl,:,iλ),
                     view(ps.Z[l],:,1:nsolve),
                     view(x, i0+1:i0+nsolve),
                     t, false)
                vl[:,iλ+1] .= conj.(vl[:,iλ])
            end
            nλ = 2
        else
            # A₁x₁ = T₁[1,1]*Z₂[:,1] = μ x₂, etc.
            il = 0
            fac = one(T)
            μ = (ps.values[1] + 0im) ^ (1/RT(p))
            μ = abs(μ)
            for l in 1:nmat
                if l == ps.schurindex
                    Tl1 = ps.T1[1,1]
                else
                    il += 1
                    Tl1 = ps.T[il][1,1]
                end
                Vs[l][:,iλ] .= fac .* ps.Z[l][:,1]
                μ == 0 || (fac *= (Tl1 / μ))
            end
            nλ = 1
        end

        sel[1:nλ] .= false
        # should we try to recover if this fails?
        ordschur!(ps, sel)
        iλ += nλ
        circshift!(sel, -nλ)
    end
    return Vs
end


# conversions
"""
     ts2fm(A::Vector{<:AbstractMatrix}, period; method = "linear") -> At::Function

Compute the function matrix corresponding to an interpolated matrix time series. 
For the given matrix time series `A`, a function matrix `A(t)` is defined as the 
mapping `A(t) = t -> Aint(t)`, where `Aint(t)` is an array of interpolation objects,  
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
   for i = 2:N
      (n1, n2) == size(A[i]) || throw(ArgumentError("all component matrices must have the same dimesnions"))
   end
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

