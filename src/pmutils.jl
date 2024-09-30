# eye(n) = Matrix{Bool}(I, n, n)
# eye(m,n) = Matrix{Bool}(I, m, n)
# eye(::Type{T}, n) where {T} = Matrix{T}(I, n, n)
eye(::Type{T}, m, n) where {T} = Matrix{T}(I, m, n)
"""
     tvstm(A, tf, t0; solver, reltol, abstol, dt) -> Φ 

Compute the state transition matrix for a linear ODE with periodic time-varying coefficients. 
For the given square periodic continuous-time matrix `A(t)`, initial time `t0` and 
final time `tf`, the state transition matrix `Φ(tf,t0)`
is computed by integrating numerically the homogeneous linear ODE 

      dΦ(t,t0)/dt = A(t)Φ(t,t0),  Φ(t0,t0) = I

on the time interval `[t0,tf]`.  
`A` may be given as a [`PeriodicFunctionMatrix`](@ref), a [`HarmonicArray`](@ref), a [`PeriodicSymbolicMatrix`](@ref)
or a  [`FourierFunctionMatrix`](@ref). 

The ODE solver to be employed can be 
specified using the keyword argument `solver` (see below), together with
the required relative accuracy `reltol` (default: `reltol = 1.e-3`), 
absolute accuracy `abstol` (default: `abstol = 1.e-7`) and/or 
the fixed step length `dt` (default: `dt = tf-t0`). 
Depending on the desired relative accuracy `reltol`, 
lower order solvers are employed for `reltol >= 1.e-4`, 
which are generally very efficient, but less accurate. If `reltol < 1.e-4`,
higher order solvers are employed able to cope with high accuracy demands. 

The following solvers from the [OrdinaryDiffEq.jl](https://github.com/SciML/OrdinaryDiffEq.jl) package can be selected:

`solver = "non-stiff"` - use a solver for non-stiff problems (`Tsit5()` or `Vern9()`);

`solver = "stiff"` - use a solver for stiff problems (`Rodas4()` or `KenCarp58()`);

`solver = "linear"` - use a special solver for linear ODEs (`MagnusGL6()`) with fixed time step `dt`;

`solver = "symplectic"` - use a symplectic Hamiltonian structure preserving solver (`IRKGL16()`);

`solver = ""` - use the default solver, which automatically detects stiff problems (`AutoTsit5(Rosenbrock23())` or `AutoVern9(Rodas5())`). 
"""
function tvstm(A::PM, tf::Real, t0::Real = 0; solver = "", reltol = 1e-3, abstol = 1e-7, dt = (tf-t0)/10) where 
         {T, PM <: Union{PeriodicFunctionMatrix{:c,T},HarmonicArray{:c,T}}} 
   n = size(A,1)
   n == size(A,2) || error("the function matrix must be square")

   isconstant(A) && ( return exp(tpmeval(A,t0)*(tf-t0)) )
   
   T1 = promote_type(typeof(t0), typeof(tf))

   # using OrdinaryDiffEq
   u0 = Matrix{T}(I,n,n)
   tspan = (T1(t0),T1(tf))
   if solver != "linear" 
      LPVODE!(du,u,p,t) = mul!(du,tpmeval(A,t),u)
      prob = ODEProblem(LPVODE!, u0, tspan)
   end
   if solver == "stiff" 
      if reltol > 1.e-4  
         # standard stiff
         sol = solve(prob, Rodas4(); reltol, abstol, save_everystep = false)
      else
         # high accuracy stiff
         sol = solve(prob, KenCarp58(); reltol, abstol, save_everystep = false)
      end
   elseif solver == "non-stiff" 
      if reltol > 1.e-4  
         # standard non-stiff
         sol = solve(prob, Tsit5(); reltol, abstol, save_everystep = false)
      else
         # high accuracy non-stiff
         sol = solve(prob, Vern9(); reltol, abstol, save_everystep = false)
      end
   elseif solver == "linear" 
      function update_func!(A,u,p,t)
         A .= p(t)
      end
      DEop = DiffEqArrayOperator(ones(T,n,n),update_func=update_func!)     
      prob = ODEProblem(DEop, u0, tspan, t-> tpmeval(A,t))
      sol = solve(prob,MagnusGL6(), dt = dt, save_everystep = false)
   elseif solver == "symplectic" 
      # high accuracy symplectic
      sol = solve(prob, IRKGaussLegendre.IRKGL16(); adaptive = true, reltol, abstol, save_everystep = false)
   else 
      if reltol > 1.e-4  
         # low accuracy automatic selection
         sol = solve(prob, AutoTsit5(Rosenbrock23()) ; reltol, abstol, save_everystep = false)
      else
         # high accuracy automatic selection
         sol = solve(prob, AutoVern9(Rodas5(),nonstifftol = 11/10); reltol, abstol, save_everystep = false)
      end
   end

   return sol(tf)     
end
""" 
     monodromy(A[, K = 1]; solver, reltol, abstol, dt) -> Ψ::PeriodicArray 

Compute the monodromy matrix for a linear ODE with periodic time-varying coefficients. 

For the given square periodic matrix `A(t)` of period `T` and subperiod `T′ = T/k`, where 
`k` is the number of subperiods,  
the monodromy matrix `Ψ = Φ(T′,0)` is computed, where `Φ(t,τ)` is the state transition matrix satisfying the homogeneous linear ODE 

    dΦ(t,τ)/dt = A(t)Φ(t,τ),  Φ(τ,τ) = I. 

`A` may be given as a [`PeriodicFunctionMatrix`](@ref), a [`HarmonicArray`](@ref), a [`PeriodicSymbolicMatrix`](@ref)
or a  [`FourierFunctionMatrix`](@ref). 

If `K > 1`, then `Ψ = Φ(T′,0)` is determined as a product of `K` matrices 
`Ψ = Ψ_K*...*Ψ_1`, where for `Δ := T′/K`, `Ψ_i = Φ(iΔ,(i-1)Δ)` is the 
state transition matrix on the time interval `[(i-1)Δ,iΔ]`. 

The state transition matrices `Φ(iΔ,(i-1)Δ)`
are computed by integrating numerically the above homogeneous linear ODE.  
The ODE solver to be employed can be 
specified using the keyword argument `solver`, together with
the required relative accuracy `reltol` (default: `reltol = 1.e-3`), 
absolute accuracy `abstol` (default: `abstol = 1.e-7`) and/or 
the fixed step length `dt` (default: `dt = min(Δ, Δ*K′/100)`) (see [`tvstm`](@ref)). 
For large values of `K`, parallel computation of factors can be alternatively performed 
by starting Julia with several execution threads. 
The number of execution threads is controlled either by using the `-t/--threads` command line argument 
or by using the `JULIA_NUM_THREADS` environment variable.  
"""
function monodromy(A::PM, K::Int = 1; solver = "non-stiff", reltol = 1e-3, abstol = 1e-7, dt = A.period/max(K,100)) where
         {T, PM <: Union{PeriodicFunctionMatrix{:c,T},HarmonicArray{:c,T}}} 
   n = size(A,1)
   n == size(A,2) || error("the periodic matrix must be square")
   nperiod = A.nperiod
   Ts = A.period/K/nperiod

   M = Array{float(T),3}(undef, n, n, K) 

   # compute the matrix exponential for K = 1 and constant matrix
   K == 1 && isconstant(A) && ( M[:,:,1] = exp(tpmeval(A,0)*Ts); return PeriodicArray(M, A.period; nperiod) )

   K >= 100 ? dt = Ts : dt = Ts*K/100/nperiod

   Threads.@threads for i = 1:K
      @inbounds M[:,:,i] = tvstm(A, i*Ts, (i-1)*Ts; solver = solver, reltol = reltol, abstol = abstol, dt = dt) 
   end
   return PeriodicArray(M,A.period; nperiod)
end
""" 
     monodromy(A) -> Ψ::PeriodicArray 

Compute the monodromy matrix for a continuous-time time series matrix. 

For the given square periodic matrix `A(t)` of period `T` and subperiod `T′ = T/k`, where 
`k` is the number of subperiods,  
the monodromy matrix `Ψ = Φ(T′,0)` is computed, where `Φ(t,τ)` is the state transition matrix satisfying the homogeneous linear ODE 

    dΦ(t,τ)/dt = A(t)Φ(t,τ),  Φ(τ,τ) = I. 

`A` is a [`PeriodicTimeSeriesMatrix`](@ref) with `K` component matrices.  
The resulting monodromy matrix `Ψ` is stored as a discrete time periodic array with `K` component matrices, of period `T` and `k` subperiods. 

`Ψ = Φ(T′,0)` is determined as a product of `K` matrices 
`Ψ = Ψ_K*...*Ψ_1`, where for `Δ := T′/K`, `Ψ_i = Φ(iΔ,(i-1)Δ)` is the 
state transition matrix on the time interval `[(i-1)Δ,iΔ]`. 
Each state transition matrix is computed as a matrix exponential  `Φ(iΔ,(i-1)Δ) = exp(A[i]*Δ)`, 
where `A[i]` is the `i`-th component matrix of the time series representation.  
For large values of `K`, parallel computation of factors can be alternatively performed 
by starting Julia with several execution threads. 
The number of execution threads is controlled either by using the `-t/--threads` command line argument 
or by using the `JULIA_NUM_THREADS` environment variable.  
"""
function monodromy(A::PM) where {T, PM <: PeriodicTimeSeriesMatrix{:c,T}} 
   n = size(A,1)
   n == size(A,2) || error("the periodic matrix must be square")
   nperiod = A.nperiod
   K = length(A)
   Ts = A.period/K/nperiod

   M = Array{float(T),3}(undef, n, n, K) 

   Threads.@threads for i = 1:K
      @inbounds M[:,:,i] = exp(A.values[i]*Ts) 
   end
   return PeriodicArray(M,A.period; nperiod)
end
function monodromy_sw(A::PM) where {T, PM <: PeriodicSwitchingMatrix{:c,T}} 
   n = size(A,1)
   n == size(A,2) || error("the periodic matrix must be square")
   nperiod = A.nperiod
   K = length(A)

   M = Array{float(T),3}(undef, n, n, K) 
   Δt = [diff(A.ts); A.period-A.ts[K]]
   Threads.@threads for i = 1:K
      @inbounds M[:,:,i] = exp(A.values[i]*Δt[i]) 
   end
   return M
end

"""
     pseig(A, K = 1; lifting = false, solver, reltol, abstol, dt) -> ev

Compute the characteristic multipliers of a continuous-time periodic matrix. 

For the given square periodic matrix `A(t)` of period `T`, 
the characteristic multipliers `ev` are the eigenvalues of 
the monodromy matrix `Ψ = Φ(T,0)`, where `Φ(t,τ)` is the state transition matrix satisfying the homogeneous linear ODE 

    dΦ(t,τ)/dt = A(t)Φ(t,τ),  Φ(τ,τ) = I. 

If `lifting = false`, `Ψ` is computed as a product of `K` state transition matrices 
`Ψ = Ψ_K*...*Ψ_1` (see [`monodromy`](@ref) with the associated keyword arguments). 
The eigenvalues are computed using the periodic Schur decomposition method of [1].
`A` may be given as a [`PeriodicFunctionMatrix`](@ref), a [`HarmonicArray`](@ref), a [`PeriodicSymbolicMatrix`](@ref)
or a  [`FourierFunctionMatrix`](@ref). 

If `lifting = true`, `Ψ` is (implicitly) expressed as `Ψ = inv(N)*M`, where `M-λN` is a regular
pencil with `N` invertible and  
the eigenvalues of `M-λN` are the same as those of the matrix product
`Ψ := Ψ_K*...*Ψ_1`. 
An efficient version of the structure exploiting fast reduction method of [2] is employed, 
which embeds the determination of transition matrices into the reduction algorithm. 
This option may occasionally lead to inaccurate results for large values of `K`. 

   _References_

[1] A. Bojanczyk, G. Golub, and P. Van Dooren, 
    The periodic Schur decomposition. Algorithms and applications, Proc. SPIE 1996.

[2] A. Varga & P. Van Dooren. Computing the zeros of periodic descriptor systems.
    Systems and Control Letters, 50:371-381, 2003.

"""
function pseig(at::PM, K::Int = 1; lifting::Bool = false, solver = "non-stiff", reltol = 1e-3, abstol = 1e-7, dt = at.period/100/at.nperiod) where 
   {T, PM <: Union{PeriodicFunctionMatrix{:c,T},HarmonicArray{:c,T}}} 
   n = size(at,1)
   n == size(at,2) || error("the periodic matrix must be square")
   nperiod = at.nperiod
   t = 0  
   Ts = at.period/K/nperiod
   if K == 1
      ev = eigvals(tvstm(at, at.period/nperiod, 0; solver, reltol, abstol, dt)) 
   else 
      if lifting  
         Z = zeros(T,n,n)
         ZI = [ Z; -I]
         si = tvstm(at, Ts, 0; solver, reltol, abstol); ti = -I
         t = Ts
         for i = 1:K-1
             tf = t+Ts
             F = qr([ ti; tvstm(at, tf, t; solver, reltol, abstol, dt) ])     
             si = F.Q'*[si; Z];  si = si[n+1:end,:]
             ti = F.Q'*ZI; ti = ti[n+1:end,:]
             t = tf
         end
         ev = -eigvals(si,ti)
         sorteigvals!(ev)
      else
         M = monodromy(at, K; solver, reltol, abstol, dt) 
         ev = K == 1 ? eigvals(view(M.M,:,:,1)) : pschur(M.M; withZ = false)[3]
         isreal(ev) && (ev = real(ev))
      end
   end
   return nperiod == 1 ? ev : ev.^nperiod
end
pseig(at::PeriodicSwitchingMatrix{:c,T}) where T = pseig(monodromy_sw(at)).^at.nperiod
pseig(at::PeriodicTimeSeriesMatrix{:c,T}) where T = pseig(monodromy(at)).^at.nperiod

"""
     ev = pseig(A::PeriodicArray; rev = true, fast = false) 

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
function pseig(A::PeriodicArray{:d,T}; fast::Bool = false) where T
   pseig(A.M; fast).^(A.nperiod)
end
function pseig(A::SwitchingPeriodicMatrix{:d,T}; fast::Bool = false) where T
   pseig(convert(PeriodicArray,A).M; fast).^(A.nperiod)
end
function pseig(A::Array{T,3}; rev::Bool = true, fast::Bool = false) where T
   n = size(A,1)
   n == size(A,2) || error("A must have equal first and second dimensions") 
   p = size(A,3)
   if fast 
      if rev 
         ev = eigvals(psreduc_reg(A)...)
      else
         imap = p:-1:1                     
         ev = eigvals(psreduc_reg(view(A,imap))...)
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
     ev = pseig(A::PeriodicMatrix[, k = 1]; rev = true, fast = false) 

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
function pseig(A::PeriodicMatrix{:d,T}, k::Int = 1; fast::Bool = false) where T
   pseig(A.M, k; fast).^(A.nperiod)
end
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
"""
     psceig(A[, K = 1]; lifting = false, solver, reltol, abstol, dt) -> ce

Compute the characteristic exponents of a continuous-time periodic matrix.

For a given square continuous-time periodic matrix `A(t)` of period `T`, 
the characteristic exponents `ce` are computed as `log.(ev)/T`, 
where  `ev` are the characteristic
multipliers (i.e., the eigenvalues of the monodromy matrix of `A(t)`).  
For available options see [`pseig(::PeriodicFunctionMatrix)`](@ref). 
`A` may be given as a [`PeriodicFunctionMatrix`](@ref), a [`HarmonicArray`](@ref), a [`PeriodicSymbolicMatrix`](@ref)
or a  [`FourierFunctionMatrix`](@ref). 
For a given square discrete-time periodic matrix `A(t)` of discrete period `N`,  
the characteristic exponents `ce` are computed as `ev.^-N`. 
"""
function psceig(at::PM, K::Int = 1; kwargs...) where
   {T, PM <: Union{PeriodicFunctionMatrix{:c,T},HarmonicArray{:c,T}}} 
   if isconstant(at)
      ce = eigvals(at(0))
   else
      ce = log.(complex(pseig(at, K; kwargs...)))/at.period
   end
   return isreal(ce) ? real(ce) : ce
end
function psceig(at::PeriodicTimeSeriesMatrix, K::Int = 1; method = "cubic", kwargs...) 
   if method == "constant"
      M = monodromy(at) 
      ev = length(at) == 1 ? eigvals(view(M.M,:,:,1)) : pschur(M.M; withZ = false)[3]
      ce = log.(complex(ev))/at.period
  else
      ce = log.(complex(pseig(convert(PeriodicFunctionMatrix,at; method), K; kwargs...)))/at.period
   end
   return isreal(ce) ? real(ce) : ce
end
function psceig(at::PeriodicSwitchingMatrix) 
   M = monodromy(at) 
   ev = length(at) == 1 ? eigvals(view(M.M,:,:,1)) : pschur(M.M; withZ = false)[3]
   ce = log.(complex(ev))/at.period
   return isreal(ce) ? real(ce) : ce
end

"""
    psceighr(Ahr::HarmonicArray[, N]; P, nperiod, shift, atol) -> ce

Compute the characteristic exponents of a continuous-time periodic matrix in _Harmonic Array_ representation. 

For a given square continuous-time periodic function matrix `Ahr(t)` of period `T` 
in a  _Harmonic Array_ representation, 
the characteristic exponents `ce` are computed as the eigenvalues of a truncated harmonic state operator `A(N)-E(N)` lying in the 
fundamental strip `-ω/2 <  Im(λ) ≤ ω/2`, where `ω = 2π/T`. If `Ahr(t)` has the harmonic components `A_0`, `A_1`, ..., `A_p`, then 
for `N ≥ p`, `P = 1` and `nperiod = 1`, the matrices `A(N)` and `E(N)` are built as


           ( A_0  A_{-1} …  A_{-p}        0    )           ( -im*ϕ_{-N}I                                 0        )
           ( A_1   A_0             ⋱           )           (     ⋮       ⋱                                        )
           (  ⋮         ⋱            ⋱         )           (               -im*ϕ_{-1}*I                           )
    A(N) = ( A_p             ⋱          A_{-p} ) , E(N) =  (                           -im*ϕ_{0}*I                )
           (        ⋱           ⋱         ⋮    )           (     ⋮                                  ⋱              )
           (  0        A_p      …         A_0  )           (     0                                   -im*ϕ_{N}I   )

with `ϕ_{i} := shift+i*ω`. If `N < p`, then a truncated _full_ block Toeplitz matrix A(N) is built using the first `N` harmonic components. 
The default value used for `N` is `N = max(10,p-1)`. 
           
Generally, for given `P ≥ 1` and  `nperiod ≥ 1`, the block Toeplitz matrix `A(N)` (and also `E(N)`) is constructed with `(2N*np+1)×(2N*np+1)` blocks,
with `np = P*nperiod`, such that each `A_i` is preceeded in its column by `np-1` zero blocks, 
each `A_{-i}` is preceeded in its row by `np-1` zero blocks and all diagonal blocks are equal to`A_0`.  

The keyword argument `atol` (default: `atol = 1.e-10`) is a tolerance on the magnitude of the trailing components of the 
associated eigenvectors used to validate their asymptotic (exponential) decay. 
Only eigenvalues satisfying this check are returned in `ce`. 

_References_

[1] J. Zhou, T. Hagiwara, and M. Araki. 
    Spectral characteristics and eigenvalues computation of the harmonic state operators in continuous-time periodic systems. 
    Systems and Control Letters, 53:141–155, 2004.
"""
function psceighr(Ahr::HarmonicArray{:c,T}, N::Int = max(10,size(Ahr.values,3)-1); P::Int = 1, nperiod::Int = Ahr.nperiod, shift::Real = 0, atol::Real = 1.e-10) where T
   n = size(Ahr,1)
   n == size(Ahr,2) || error("the periodic matrix must be square") 
   (N == 0 || isconstant(Ahr)) && (return eigvals(real(Ahr.values[:,:,1])))
   ev, V = eigen!(hr2btupd(Ahr, N; P, shift, nperiod));
   ind = sortperm(imag(ev),by=abs); 
   ωhp2 = pi/P/Ahr.period/Ahr.nperiod*nperiod
   ne = count(abs.(imag(ev[ind[1:min(4*n,length(ind))]])) .<=  ωhp2*(1+sqrt(eps(T))))
   ev = ev[ind[1:ne]]
   # return only validated eigenvalues
   σ = Complex{T}[]
   for i = 1:ne
       norm(V[end-3:end,ind[i]]) < atol  && push!(σ,ev[i])
   end
   nv = length(σ)
   nv < n && @warn "number of eigenvalues is less than the order of matrix, try again with increased number of harmonics"
   ce = nv > n ? σ[sortperm(imag(σ),rev=true)][1:n] : σ[1:nv]
   return isreal(ce) ? real(ce) : ce
end
"""
    psceig(A::AbstractPeriodicArray[, k]; kwargs...) -> ce

Compute the characteristic exponents of a cyclic matrix product of `p` matrices.

The characteristic exponents of a product of `p` matrices are computed as the `p`th roots of the 
characteristic multipliers. These are computed as the eigenvalues of the square 
cyclic product of `p` matrices `A(k-1)...*A(2)*A(1)*A(p)...*A(k)`, if `rev = true` (default) or 
`A(k)*A(k+1)*...A(p)*A(1)...A(k-1)` if `rev = false`, without evaluating the product. 
The argument `k` specifies the starting index (default: `k = 1`). 
The matrices `A(1)`, `...`, `A(p)` are contained in the `p`-vector of matrices `A` 
such that the `i`-th matrix  `A(i)`, of dimensions `m(i)×n(i)`, is contained in `A[i]`.
The keyword arguments `kwargs` are those of  [`pseig(::PeriodicMatrix)`](@ref).  

_Note:_ The first `nmin` components of `ce` contains the _core characteristic exponents_ of the appropriate matrix product,
where `nmin` is the minimum row dimensions of matrices `A[i]`, for `i = 1, ..., p`, 
while the last components of `ce` are zero. 
"""
function psceig(at::AbstractPeriodicArray{:d,T}, k::Int = 1; kwargs...) where T
   if isconstant(at)
      ce = eigvals(at(1))
   else
      ce = (complex(pseig(convert(PeriodicMatrix,at), k; kwargs...))).^(1/at.dperiod/at.nperiod) 
   end
   return isreal(ce) ? real(ce) : ce
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
     ts2pfm(At::PeriodicTimeSeriesMatrix; method = "linear") -> A::PeriodicFunctionMatrix

Compute the periodic function matrix corresponding to an interpolated periodic time series matrix. 
For the given periodic time series matrix `At`, a periodic function matrix `A(t)` is defined as the 
mapping `A(t) = t -> etpf(t)`, where `etpf(t)` is a periodic interpolation/extrapolation object,  
as provided in the [`Interpolations.jl`](https://github.com/JuliaMath/Interpolations.jl)  package. 
The keyword parameter `method` specifies the interpolation/extrapolation method to be used as follows:

`method = "constant"` - use periodic B-splines of degree 0 (periodic constant interpolation);

`method = "linear"` - use periodic B-splines of degree 1 (periodic linear interpolation) (default);

`method = "quadratic"` - use periodic B-splines of degree 2 (periodic quadratic interpolation); 

`method = "cubic"` - use periodic B-splines of degree 3 (periodic cubic interpolation). 
"""
function ts2pfm(A::PeriodicTimeSeriesMatrix; method = "linear")
   N = length(A.values)
   N == 0 && error("empty time array not supported")
   N == 1 && (return PeriodicFunctionMatrix(t -> A.values[1], A.period; nperiod = A.nperiod, isconst = true))
   #ts = (0:N-1)*(A.period/N)
   d = A.period/N/A.nperiod
   ts = (0:N-1)*d
   n1, n2 = size(A.values[1])
   intparray = Array{Any,2}(undef, n1, n2)
   if method == "constant"     
      # use simple function evaluation
      return PeriodicFunctionMatrix(t -> tpmeval(A,t), A.period; nperiod = A.nperiod, isconst = isconstant(A))
      # ts = ts.+d/2  # use center points
      # [intparray[i,j] = scale(Interpolations.extrapolate(interpolate(getindex.(A.values,i,j), BSpline(Constant(Periodic(OnCell())))), Periodic()), ts) for i in 1:n1, j in 1:n2]
   elseif method == "linear" || N == 2
      [intparray[i,j] = scale(Interpolations.extrapolate(interpolate(getindex.(A.values,i,j), BSpline(Linear(Periodic(OnCell())))), Periodic()), ts) for i in 1:n1, j in 1:n2]
   elseif method == "quadratic" || N == 3     
      [intparray[i,j] = scale(Interpolations.extrapolate(interpolate(getindex.(A.values,i,j), BSpline(Quadratic(Periodic(OnCell())))), Periodic()), ts) for i in 1:n1, j in 1:n2]
   elseif method == "cubic"     
      [intparray[i,j] = scale(Interpolations.extrapolate(interpolate(getindex.(A.values,i,j), BSpline(Cubic(Periodic(OnCell())))), Periodic()), ts) for i in 1:n1, j in 1:n2]
   else
      error("no such option method = $method")
   end
   return PeriodicFunctionMatrix(t -> [intparray[i,j](t) for i in 1:n1, j in 1:n2 ], A.period; nperiod = A.nperiod, isconst = isconstant(A))
end
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
      #[intparray[i,j] = scale(Interpolations.extrapolate(interpolate(getindex.(A,i,j), BSpline(Cubic(OnGrid()))), Interpolations.Flat()), ts) for i in 1:n1, j in 1:n2]
      #[intparray[i,j] = scale(interpolate(getindex.(A,i,j), BSpline(Cubic(Interpolations.Line(OnGrid())))), ts) for i in 1:n1, j in 1:n2]
      [intparray[i,j] = scale(interpolate(getindex.(A,i,j), BSpline(Cubic())), ts) for i in 1:n1, j in 1:n2]
   else
      error("no such option method = $method")
   end
   return t -> [intparray[i,j](t) for i in 1:n1, j in 1:n2 ]
end

# function tsw2pfm(A::PeriodicSwitchingMatrix; method = "linear")
#    N = length(A.values)
#    N == 0 && error("empty time array not supported")
#    N == 1 && (return PeriodicFunctionMatrix(t -> A.values[1], A.period; nperiod = A.nperiod, isconst = true))
#    #ts = (0:N-1)*(A.period/N)
#    # d = A.period/N/A.nperiod
#    # ts = (0:N-1)*d
#    ts = A.ts
#    n1, n2 = size(A.values[1])
#    intparray = Array{Any,2}(undef, n1, n2)
#    if method == "constant"     
#       # use simple function evaluation
#       return PeriodicFunctionMatrix(t -> tpmeval(A,t), A.period; nperiod = A.nperiod, isconst = isconstant(A))
#       # ts = ts.+d/2  # use center points
#       # [intparray[i,j] = scale(Interpolations.extrapolate(interpolate(getindex.(A.values,i,j), BSpline(Constant(Periodic(OnCell())))), Periodic()), ts) for i in 1:n1, j in 1:n2]
#    elseif method == "linear" || N == 2
#       [intparray[i,j] = linear_interpolation(ts,getindex.(A.values,i,j); extrapolation_bc=Periodic()) for i in 1:n1, j in 1:n2]
#    else
#       error("no such option method = $method")
#    end
#    return PeriodicFunctionMatrix(t -> [intparray[i,j](t) for i in 1:n1, j in 1:n2 ], A.period; nperiod = A.nperiod, isconst = isconstant(A))
# end

"""
     ts2hr(A::PeriodicTimeSeriesMatrix; atol = 0, rtol, n, squeeze = true) -> Ahr::HarmonicArray

Compute the harmonic (Fourier) approximation of a periodic matrix specified by a time series matrix. 
The periodic matrix `A(t)` is specified as a continuous-time periodic time series matrix `A`, 
with `m` matrices contained in the vector of matrices `A.values`, where `A.values[k]` 
is the value of `A(t)` at time moment `(k-1)T/m`, with `T = A.period` being the period. 
The resulting harmonic approximation `Ahr(t)` of `A(t)` has the form

                     p
     Ahr(t) = A_0 +  ∑ ( Ac_i*cos(i*t*2*π/T)+As_i*sin(i*2*π*t/T) ) 
                    i=1 

where `A_0` is the constant term (the mean value), `Ac_i` and `As_i` are the  
coefficient matrices of the `i`-th cosinus and sinus terms, respectively. 
The order of the approximation `p` is determined using the maximum order specified by `n` 
(default: `n = (m-1)/2`) and the absolute and relative tolerances `atol` and `rtol`, as follows:
`p` is the minimum between `n`, `(m-1)/2` and the maximum index `k` 
such that `Ac_k` and/or `As_k` are nonzero.
The tolerance used to assess nonzero elements is `tol = max(atol, rtol*maxnorm)`, where 
`maxnorm` is the maximum norm of the matrices contained in `A.values`. The default values of tolerances
are `atol = 0` and `rtol = 10*p*ϵ`, where `ϵ` is the working machine precision.

The resulting harmonic approximation `Ahr(t)` is returned in the harmonic array object `Ahr` 
(see [`HarmonicArray`](@ref)). 
"""
function ts2hr(A::PeriodicTimeSeriesMatrix{:c,T}; atol::Real = 0, rtol::Real = 0, n::Union{Int,Missing} = missing, squeeze::Bool = true) where  {T}
   
   M = length(A.values)
   n1, n2 = size(A.values[1])
   
   if ismissing(n)
       n = div(M-1,2)
       ncur = 0
   else
       n = min( n, div(M-1,2) ) 
       ncur = n
   end
   n = max(n,0)
   
   AHR = zeros(ComplexF64, n1, n2, n+1)
   T1 = promote_type(Float64,T)
   tol = iszero(atol) ? (iszero(rtol) ? 10*M*eps(T1)*maximum(norm.(A.values)) : rtol*maximum(norm.(A.values)) ) : atol
   i1 = 1:n+1   
   for i = 1:n1
       for j = 1:n2
           temp = T1.(getindex.(A.values, i, j))
           i == 1 && j == 1 && (global rfftop = plan_rfft(temp; flags=FFTW.ESTIMATE, timelimit=Inf))
           tt = conj(2/M*(rfftop*temp)) 
           tt[1] = real(tt[1]/2)
           tt1 = view(tt,i1)
           indr = i1[abs.(real(tt1)) .> tol]
           nr = length(indr); 
           nr > 0 && (nr = indr[end])
           indi = i1[abs.(imag(tt1)) .> tol]
           ni = length(indi); 
           ni > 0 && (ni = indi[end])
           ncur = max(ncur, nr, ni)        
           AHR[i,j,indr] = real(tt[indr])
           AHR[i,j,indi] += im*imag(tt[indi])
       end
   end
   nperiod0 = A.nperiod
   nperiod = 1
   if ncur > 2 && squeeze
      nh = ncur-1
      s = falses(nh)
      for i = 1:nh
          s[i] = any(abs.(view(AHR,:,:,i+1)) .> tol)
      end  
      t = Primes.factor(Vector,nh)
      s1 = copy(s)
     for i = 1:length(t)
          stry = true
          for j1 = 1:t[i]:nh
              stry = stry & all(view(s1,j1:j1+t[i]-2) .== false) 
              stry || break
          end
          if stry 
             nperiod = nperiod*t[i]
             s1 = s1[t[i]:t[i]:nh]
             nh = div(nh,t[i])
          end
      end 
      return HarmonicArray(AHR[:,:,[1;nperiod+1:nperiod:ncur]],A.period;nperiod = nperiod*nperiod0)
   else
      return HarmonicArray(AHR[:,:,1:max(1,ncur)],A.period;nperiod = nperiod0)
   end       

end
"""
     pfm2hr(A::PeriodicFunctionMatrix; nsample, NyquistFreq) -> Ahr::HarmonicArray

Convert a periodic function matrix into a harmonic array representation. 
If `A(t)` is a periodic function matrix of period `T`, then the harmonic array representation
`Ahr` is determined by applying the fast Fourier transform to the sampled arrays `A(iΔ)`, `i = 0, ..., k`,
where `Δ = T/k` is the sampling period and `k` is the number of samples specified by the keyword argument
`nsample = k` (default: `k = 128`). If the Nyquist frequency `f` is specified via the keyword argument
`NyquistFreq = f`, then `k` is chosen `k = 2*f*T` to avoid signal aliasing.     
"""
function pfm2hr(A::PeriodicFunctionMatrix; nsample::Int = 128, NyquistFreq::Union{Real,Missing} = missing)   
   isconstant(A) && (return HarmonicArray(A.f(0),A.period; nperiod = A.nperiod))
   nsample > 0 || ArgumentError("nsample must be positive, got $nsaple")
   ns = ismissing(NyquistFreq) ? nsample : Int(floor(2*abs(NyquistFreq)*A.period/A.nperiod))+1
   Δ = A.period/ns/A.nperiod
   ts = (0:ns-1)*Δ
   return ts2hr(PeriodicTimeSeriesMatrix(A.f.(ts), A.period; nperiod = A.nperiod),squeeze = true)
end
"""
     hrchop(Ahr::HarmonicArray; tol) -> Ahrtrunc::HarmonicArray

Remove the trailing terms of a harmonic representation by deleting those whose norms are below a certain tolerance. 
"""
function hrchop(ahr::HarmonicArray; tol::Real = sqrt(eps()) ) 
   nrm = norm(ahr.values[:,:,1],1)
   n3 = size(ahr.values,3)
   for i = 2:n3
       nrm = max(nrm,norm(ahr.values[:,:,i],1))
   end
   itrunc = 1
   atol = tol*nrm
   for i = n3:-1:1
       if norm(ahr.values[:,:,i],1) > atol 
          itrunc = i
          break
       end
   end
   return HarmonicArray(ahr.values[:,:,1:itrunc], ahr.period; nperiod = ahr.nperiod)
end
"""
     hrtrunc(Ahr::HarmonicArray, n) -> Ahrtrunc::HarmonicArray

Truncate a harmonic representation by deleting the trailing terms whose indices exceed certain number `n` of harmonics. 

"""
function hrtrunc(ahr::HarmonicArray, n::Int = size(ahr.values,3)-1) 
   return HarmonicArray(ahr.values[:,:,1:max(1,min(n+1,size(ahr.values,3)))], ahr.period; nperiod = ahr.nperiod)
end

"""
     pm2pa(At::PeriodicMatrix) -> A::PeriodicArray

Convert a discrete-time periodic matrix object into a discrete-time periodic array object.

The discrete-time periodic matrix object `At` contains a  
`p`-vector `At` of real matrices `At[i]`, `i = 1,..., p`, 
the associated time period `T` and the number of subperiods `k`. The resulting discrete-time periodic array object
`A` of period `T` and number of subperiods `k` 
is built from a `m×n×p` real array `A`, such that `A[:,:,i]` 
contains `At[i]`, for `i = 1,..., p`. 
For non-constant dimensions, the resulting `m` and `n` are the maximum row and column dimensions, respectively,
and the resulting component matrices `A[:,:,i]` contain `At[i]`, appropriately padded with zero entries. 
"""
function pm2pa(A::PeriodicMatrix{:d,T}) where T
   N = length(A)
   m, n = size(A)
   N == 0 && PeriodicArray(Array{T,3}(undef,0,0,0),A.period)
   if any(m .!= m[1]) || any(m .!= m[1]) 
      #@warn "Non-constant dimensions: the resulting component matrices padded with zeros"
      t = zeros(T,maximum(m),maximum(n),N)
      [copyto!(view(t,1:m[i],1:n[i],i),A.M[i]) for i in 1:N]
      PeriodicArray{:d,T}(t,A.period,A.nperiod)
   else
      t = zeros(T,m[1],n[1],N)
      [copyto!(view(t,:,:,i),A.M[i]) for i in 1:N]
      PeriodicArray{:d,T}(t,A.period,A.nperiod)
   end
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

# time response evaluations

"""
     tvmeval(At::PeriodicTimeSeriesMatrix, t; method = "linear") -> A::Vector{Matrix}

Evaluate the time response of a periodic time series matrix.

For the periodic time series matrix `At` and the vector of time values `t`, 
an interpolation/extrapolation based approximation  
`A[i]` is evaluated for each time value `t[i]`. The keyword parameter `method` specifies the
interpolation/extrapolation method to be used for periodic data. 
The following interpolation methods from the [`Interpolations.jl`](https://github.com/JuliaMath/Interpolations.jl) 
package can be selected: 

`method = "constant"` - use periodic B-splines of degree 0; 

`method = "linear"` - use periodic B-splines of degree 1 (periodic linear interpolation) (default);

`method = "quadratic"` - use periodic B-splines of degree 2 (periodic quadratic interpolation); 

`method = "cubic"` - use periodic B-splines of degree 3 (periodic cubic interpolation).
"""
function tvmeval(A::PeriodicTimeSeriesMatrix{:c,T}, t::Union{Real,Vector{<:Real}}; method = "linear") where T
   N = length(A.values)
   N == 0 && error("empty time array not supported")
   isa(t,Real) ? te = [t] : te = t
   N == 1 && (return [A.values[1] for i in te])
   nt = length(te)
   nt == 0 && (return zeros(T,size(A.values,1),size(A.values,2),0))
   dt = A.period/N
   ts = (0:N-1)*dt
   n1, n2 = size(A.values[1])
   intparray = Array{Any,2}(undef,n1, n2)
   if method == "linear"
      [intparray[i,j] = scale(Interpolations.extrapolate(interpolate(getindex.(A.values,i,j), BSpline(Linear(Periodic(OnCell())))), Periodic()), ts) for i in 1:n1, j in 1:n2]
   elseif method == "cubic"      
      [intparray[i,j] = scale(Interpolations.extrapolate(interpolate(getindex.(A.values,i,j), BSpline(Cubic(Periodic(OnCell())))), Periodic()), ts) for i in 1:n1, j in 1:n2]
   elseif method == "quadratic"      
      [intparray[i,j] = scale(Interpolations.extrapolate(interpolate(getindex.(A.values,i,j), BSpline(Quadratic(Periodic(OnCell())))), Periodic()), ts) for i in 1:n1, j in 1:n2]
   elseif method == "constant"      
      [intparray[i,j] = scale(Interpolations.extrapolate(interpolate(getindex.(A.values,i,j), BSpline(Constant(Periodic(OnCell())))), Periodic()), ts) for i in 1:n1, j in 1:n2]
   else
      error("no such option method = $method")
   end
   return [[intparray[i,j].(te[k]) for i in 1:n1, j in 1:n2 ] for k = 1:nt ]
end
function tvmeval(A::PeriodicSwitchingMatrix{:c,T}, t::Union{Real,Vector{<:Real}}) where T
   [tpmeval(A, x) for x in t]
end

"""
     hreval(Ahr::HarmonicArray, t; ntrunc, exact = true) -> A::Matrix

Evaluate the harmonic array `Ahr` representing a continuous-time 
time periodic matrix `A(t)` for a time value `t`. 
For a real value `t`, if `exact = true (default)` an exact evaluation is computed, while for `exact = false`, 
a linear interpolation based approximation is computed (potentially more accurate in intersample points).
The keyword argument `ntrunc` specifies the number of harmonics to be used for the evaluation 
(default: maximum possible number). 
"""
function hreval(ahr::HarmonicArray{:c,T}, t::Real; exact::Bool = true, ntrunc::Int = max(size(ahr.values,3)-1,0)) where {T}
      (ntrunc < 0 || ntrunc >= size(ahr.values,3)) && error("ntrunc out of allowed range")

   n = ntrunc   
   T1 = float(promote_type(T,typeof(t)))
   ts = mod(T1(t),T1(ahr.period))*2*pi*ahr.nperiod/ahr.period
   
   # determine interpolation coefficients
   ht = ones(T1,n);
   if !exact
      # use linear interpolation
      for i = 2:n
           x = pi*(i-1)/n
           ht[i] = (sin(x)/x)^2
      end
   end
   a = T == T1 ? real.(ahr.values[:,:,1]) : T1.(real.(ahr.values[:,:,1]))
   for i = 1:n
       ta = view(ahr.values,:,:,i+1)
       a .+= T1.(real.(ta)).*(cos(i*ts)*ht[i]) .+ T1.(imag.(ta)) .* ((sin(i*ts)*ht[i]))
   end
   return a
end   
"""
     tvmeval(Ahr::HarmonicArray, t; ntrunc, exact = true) -> A::Vector{Matrix}

Evaluate the time response of a harmonic array.

For the harmonic array `Ahr` representing representing a continuous-time 
time periodic matrix `A(t)` and the vector of time values `t`, 
`A[i] = A(t[i])` is computed for each time value `t[i]`. 
If `exact = true (default)` an exact evaluation is computed, while for `exact = false`, 
a linear interpolation based approximation is computed 
(potentially more accurate in intersample points).
The keyword argument `ntrunc` specifies the number of harmonics to be used for evaluation 
(default: maximum possible number of harmonics). 
"""
function tvmeval(ahr::HarmonicArray{:c,T}, t::Union{Real,Vector{<:Real}}; ntrunc::Int = size(ahr.values,3), 
                exact::Bool = true) where {T}
       
   n = min(size(ahr.values,3)-1,ntrunc);
   
   isa(t,Real) ? te = [t] : te = t
   nt = length(te)
   period = ahr.period
   
   tscal = 2*pi*ahr.nperiod/period
   # determine interpolation coefficients
   ht = ones(Float64,n);
   if !exact
      # use linear interpolation
      for i = 2:n
           x = pi*(i-1)/n
           ht[i] = (sin(x)/x)^2
      end
   end
   T1 = float(T)
   A = similar(Vector{Matrix{T1}}, nt)
   for j = 1:nt
       A[j] = real(ahr.values[:,:,1])
       tsj = mod(te[j],period)*tscal
       for i = 1:n
           ta = view(ahr.values,:,:,i+1)
           A[j] .+= real.(ta).*(cos(i*tsj)*ht[i]) .+ imag.(ta) .* ((sin(i*tsj)*ht[i]))
       end
   end
   return A
end   
"""
     tvmeval(A, t) -> Aval::Vector{Matrix}

Evaluate the time response of a periodic matrix.

For the periodic matrix `A(t)` and the vector of time values `t`, 
the vector `Aval` of time values is computed such that 
`Aval[i] = A(t[i])` for each time value `t[i]`. 
"""
function tvmeval(A::PeriodicFunctionMatrix, t::Union{Real,Vector{<:Real}} )
   return tpmeval.(Ref(A),t)
end
function tpmeval(A::PeriodicFunctionMatrix, t::Real )
"""
     tpmeval(A, tval) -> Aval::Matrix

Evaluate the time value of a periodic matrix.

For the periodic matrix `A(t)` and the time value `tval`, `Aval = A(tval)` is evaluated for `t = tval`. 
"""
   return reshape((A.f).(mod(t,A.period/A.nperiod)),size(A,1),size(A,2))
end

tpmeval(A::HarmonicArray, t::Real ) = hreval(A, t; exact = true) 
function tpmeval(A::PeriodicTimeSeriesMatrix, t::Real )
   tsub = A.period/A.nperiod
   ns = length(A.values)
   Δ = tsub/ns
   ind = floor(Int,mod(t,tsub)/Δ)+1
   ind <= ns || (ind = 1) 
   return A.values[ind]
end
function tpmeval(A::PeriodicSwitchingMatrix, t::Real )
   ind = findfirst(A.ts .> mod(t,A.period/A.nperiod)*(1+10*eps()))
   isnothing(ind) ? ind = length(A) : ind -= 1
   return A.values[ind]
end
# function tpmeval(A::SwitchingPeriodicMatrix, t::Real )
#    ts = [0; Int.(round.(A.ns[1:end-1]))]*A.Ts
#    ind = findfirst(ts .> mod(t,A.period/A.nperiod))
#    isnothing(ind) ? ind = length(A.M) : ind -= 1
#    return A.M[ind]
# end
function tpmeval(A::SwitchingPeriodicMatrix, t::Real )
   nv = length(A.ns)
   k = Int(floor(mod(t,A.period/A.nperiod) / A.Ts)) 
   ind = findfirst(view(A.ns,1:nv-1) .> mod(k,A.dperiod))
   isnothing(ind) && (ind = nv)
   return A.M[ind]
end
function tpmeval(A::PeriodicMatrix, t::Real )
   nv = length(A.M)
   k = Int(floor(mod(t,A.period/A.nperiod) / A.Ts)) 
   ind = mod(k+1,A.dperiod) 
   ind == 0 && (ind = nv)
   return A.M[ind]
end
function tpmeval(A::PeriodicArray, t::Real )
   nv = size(A.M,3)
   k = Int(floor(mod(t,A.period/A.nperiod) / A.Ts)) 
   ind = mod(k+1,A.dperiod) 
   ind == 0 && (ind = nv)
   return A.M[:,:,ind]
end
function kpmeval(A::SwitchingPeriodicMatrix, k::Int)
   ind = findfirst(A.ns .>= mod(k-1,A.dperiod)+1)
   return A.M[ind]
end
function kpmeval(A::SwitchingPeriodicArray, k::Int)
   ind = findfirst(A.ns .>= mod(k-1,A.dperiod)+1)
   return A.M[:,:,ind]
end
(F::PeriodicFunctionMatrix)(t) = tpmeval.(Ref(F), t)
(F::HarmonicArray)(t) = hreval(F, t; exact = true) 
(F::PeriodicTimeSeriesMatrix)(t) = tpmeval(F, t) 
(F::PeriodicSwitchingMatrix)(t) = tpmeval(F, t) 
(F::SwitchingPeriodicMatrix)(t) = tpmeval(F, t) 
(F::PeriodicMatrix)(t) = tpmeval(F, t) 
(F::PeriodicArray)(t) = tpmeval(F, t) 
   
"""
    pmaverage(A; rtol = sqrt(eps())) -> Am 

Compute for the continuous-time periodic matrix `A(t)` 
the corresponding time averaged matrix `Am` over one period. 
The Gauss-Konrod quadratue method is employed for numerical integration
using a relative accuracy tolerance specified by `rtol`. 
"""
function pmaverage(A::PM; rtol = sqrt(eps())) where {PM <: PeriodicFunctionMatrix} 
   tsub = A.period/A.nperiod
   tt, = quadgk(A.f, 0., tsub; rtol)
   return tt/tsub
end
function pmaverage(A::PM; rtol = sqrt(eps())) where {PM <: Union{PeriodicSwitchingMatrix,PeriodicTimeSeriesMatrix}} 
   return pmaverage(convert(PeriodicFunctionMatrix,A); rtol)
end
pmaverage(A::HarmonicArray; rtol = missing) = real(A.values[:,:,1])
function getpm(A::PeriodicMatrix, k, dperiod::Union{Int,Missing} = missing)
   i = ismissing(dperiod) ? mod(k-1,A.dperiod)+1 : mod(k-1,dperiod)+1
   return A.M[i]
   #return view(A.M,i)
end
function getpm(A::SwitchingPeriodicMatrix, k, dperiod::Union{Int,Missing} = missing)
   # i = ismissing(dperiod) ? mod(k-1,A.dperiod)+1 : mod(k-1,dperiod)+1
   return A[k]
   #return view(A.M,i)
end
function getpm(A::PeriodicArray, k, dperiod::Union{Int,Missing} = missing)
   i = ismissing(dperiod) ? mod(k-1,A.dperiod)+1 : mod(k-1,dperiod)+1
   return A.M[:,:,i]
   #return view(A.M,:,:,i)
end
function copypm!(Dest::AbstractMatrix{T}, A::PeriodicMatrix{:d,T}, k, dperiod::Union{Int,Missing} = missing) where {T}
   i = ismissing(dperiod) ? mod(k-1,A.dperiod)+1 : mod(k-1,dperiod)+1
   return copyto!(Dest,view(A.M[i],:,:))
   #return copyto!(Dest,A.M[i])
   #return view(A.M,i)
end
function copypm!(Dest::AbstractMatrix, A::PeriodicArray, k, dperiod::Union{Int,Missing} = missing)
   i = ismissing(dperiod) ? mod(k-1,A.dperiod)+1 : mod(k-1,dperiod)+1
   return copyto!(Dest,view(A.M,:,:,i))
   #return view(A.M,:,:,i)
end

"""
     hr2bt(Ahr::HarmonicArray, N; P, nperiod]) -> Abt::Matrix 

Build the block Toeplitz matrix of a harmonic (Fourier) array representation of a periodic matrix. 

The harmonic representation object `Ahr` of period `T` of a periodic matrix `Ahr(t)` 
of subperiod `T′ = T/k` is in the form

                     p
     Ahr(t) = A_0 +  ∑ ( Ac_i*cos(i*t*2*π/T′)+As_i*sin(i*2*π*t/T′) ) ,
                    i=1 

where `k ≥ 1` is the number of subperiods. `Ahr(t)` can be equivalently expressed in the Fourier series
representation form

                p
     Ahr(t) =   ∑ A_i*exp(im*i*ωh*t) ,
               i=-p

where `ωh = 2π/T′`, `A_i = (Ac_i-im*As_i)/2` and  `A_{-i} = (Ac_i+im*As_i)/2`. 
`N` is the number of selected harmonic components (or Fourier modes) used for approximation. 
The keyword parameter `P` is the number of full periods to be considered (default: `P = 1`) and `nperiod` is the
number of subperiods to be considered, such that `1 ≤ nperiod ≤ k` (default: `nperiod = k`). 

For a given number `N ≥ p`, if the number of period is `P = 1` and the number of subperiods is `nperiod = 1`, 
then the _banded_ block Toeplitz matrix `Abt` with `(2N+1)×(2N+1)` blocks is constructed
           
           ( A_0  A_{-1} …  A_{-p}        0    )
           ( A_1   A_0             ⋱           )
           (  ⋮         ⋱            ⋱         )
     Abt = ( A_p             ⋱          A_{-p} )
           (        ⋱           ⋱         ⋮    )
           (  0        A_p      …         A_0  )

If `N < p`, then a truncated _full_ block Toeplitz matrix is built using the first `N` harmonic components. 

Generally, for given `P ≥ 1` and  `nperiod ≥ 1`, the block Toeplitz matrix `Abt` is constructed with `(2N*np+1)×(2N*np+1)` blocks,
with `np = P*nperiod`, such that each `A_i` is preceeded in its column by `np-1` zero blocks, 
each `A_{-i}` is preceeded in its row by `np-1` zero blocks and all diagonal blocks are equal to`A_0`.   
"""
function hr2bt(A::HarmonicArray, N::Int; P::Int = 1, nperiod::Int = A.nperiod, nh::Int = size(A.values,3)-1)
    N > 0 || error("the number of harmonic components must be nonnegative, got $N")
    #(nperiod < 1 || nperiod > A.nperiod) && error("number of subperiods must be between 1 and $(A.nperiod), got $nperiod")
    nperiod < 1  && throw(ArgumentError("number of subperiods must be between 1 and $(A.nperiod), got $nperiod"))
    P < 1 && throw(ArgumentError("number of periods must be at least 1, got $P"))
    nh < 0 && throw(ArgumentError("number of harmonics must be non-negative, got $nh"))
    T = promote_type(Float64,eltype(A))
    # p = size(A.values,3)-1
    # p < 0 && (return zeros(complex(T),0,0))
    p = min(nh,size(A.values,3)-1)
    m, n = size(A)
    np = P*nperiod
    nb = 2*N*np+1
    BT = zeros(complex(T),nb*m,nb*n)
    i1, i2, j1, j2 = 1, m, 1, n
    BT[i1:i2,j1:j2] = A.values[:,:,1] # the diagonal A0
    ki1 = m*np+1
    ki2 = ki1+m-1
    kj1 = n*np+1
    kj2 = kj1+n-1
    minpN = min(p,N)
    for k = 1:minpN
        BT[i1:i2,kj1:kj2] = A.values[:,:,k+1]/2           # this is A_{-k} := (Ac_k+im*As_k)/2 
        BT[ki1:ki2,j1:j2] = conj(view(BT,i1:i2,kj1:kj2))  # this is A_k    := (Ac_k-im*As_k)/2 
        ki1 += m*np
        ki2 += m*np
        kj1 += n*np
        kj2 += n*np
    end
    i1 = 1
    j1 = n+1
    ik1 = m+1
    ik2 = 2m
    jk1 = n+1
    jk2 = 2n
    for i = 1:nb-1
        i1 += m
        i2 = min(i1+m*(minpN*np+1)-1,nb*m)
        i1 <= i2 || continue
        j1 += n
        j2 = min(j1+n*np-1,nb*n)
        BT[i1:i2,jk1:jk2] = BT[1:min(i2-i1+1,nb*m),1:n]
        BT[ik1:ik2,j1:j2] = BT[1:m,n+1:n+min(j2-j1+1,nb*n)]
        ik1 += m
        ik2 += m
        jk1 += n
        jk2 += n
    end
    return BT
end
"""
     hr2btupd(Ahr::HarmonicArray, N; P, nperiod, shift]) -> Ab::Matrix 

Build the updated block Toeplitz matrix of a harmonic (Fourier) array representation of a periodic matrix. 

If `BT` is the block Toeplitz matrix of the harmonic array representation of the `n × n` periodic matrix `Ahr` of subperiod `T′` 
(see [`HarmonicArray`](@ref)) as constructed with [`hr2bt`](@ref), then the updated matrix Ab = BT-NT is constructed, 
with `NT` a block-diagonal matrix with `n × n` diagonal blocks.
The `k`-th diagonal block of `NT` is the diagonal matrix `im*(μ + k*ωh)*I`, where `μ` is a shift specified via 
the keyword parameter `shift = μ` (default: `μ = 0`)  and `ωh` is the base frequency computed as `ωh = 2π*nperiod/(P*T′)`. 
The value of shift must satisfy `0 ≤ μ ≤ ωh/2`. 
"""
function hr2btupd(A::HarmonicArray, N::Int; P::Int = 1, nperiod::Int = A.nperiod, shift::Real = 0)
    n = size(A,1)
    n == size(A,2) || error("the periodic matrix must be square") 
    BT = hr2bt(A, N; P, nperiod)
    np = P*nperiod
    nb = 2*N*np+1
    ωh = 2*pi/P/A.period/A.nperiod*nperiod
    (shift >= 0 && shift <= ωh/2) 
    #Ej0 = fill(eltype(BT)(shift + im*ωh),n)
    Ej0 = similar(Vector{eltype(BT)},n)
    k = -N*np
    for i = 1:nb
        ind = (i-1)*n+1:i*n
        BT0 = view(BT,ind,ind)
        #BT0[diagind(BT0)] = diag(BT0)-k*Ej0;
        BT0[diagind(BT0)] = diag(BT0)-fill!(Ej0,(im*(shift + k*ωh)))
        k += 1
    end
    return BT
end
