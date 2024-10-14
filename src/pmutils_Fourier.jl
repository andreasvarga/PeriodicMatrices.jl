"""
     tvstm(A::FourierFunctionMatrix, tf, t0; solver, reltol, abstol, dt) -> Φ 

Compute the state transition matrix for a linear ODE with periodic time-varying coefficients. 
For the given periodic square matrix `A(t)`, initial time `t0` and 
final time `tf`, the state transition matrix `Φ(tf,t0)`
is computed by integrating numerically the homogeneous linear ODE 

      dΦ(t,t0)/dt = A(t)Φ(t,t0),  Φ(t0,t0) = I

on the time interval `[t0,tf]`. `A(t)` has the type `FourierFunctionMatrix`. 

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
         {T, PM <: FourierFunctionMatrix{:c,T}} 
   n = size(A,1)
   n == size(A,2) || error("the function matrix must be square")

   isconstant(A) && ( return exp(tpmeval(A,t0)*(tf-t0)) )
   
   T1 = promote_type(typeof(t0), typeof(tf))

   # using OrdinaryDiffEq
   u0 = Matrix{eltype(A) == Num ? Float64 : T}(I,n,n)
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
      #prob = ODEProblem(DEop, u0, tspan, A.f)
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
     monodromy(A::FourierFunctionMatrix[, K = 1]; solver, reltol, abstol, dt) -> Ψ::PeriodicArray 

Compute the monodromy matrix for a linear ODE with periodic time-varying coefficients. 

For the given square periodic matrix `A(t)` of period `T` and subperiod `T′ = T/k`, where 
`k` is the number of subperiods,  
the monodromy matrix `Ψ = Φ(T′,0)` is computed, where `Φ(t,τ)` is the state transition matrix satisfying the homogeneous linear ODE 

    dΦ(t,τ)/dt = A(t)Φ(t,τ),  Φ(τ,τ) = I. 

`A(t)` has the type `FourierFunctionMatrix`. 

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
         {T, PM <: FourierFunctionMatrix{:c,T}}
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
     pseig(A::FourierFunctionMatrix[, K = 1]; lifting = false, solver, reltol, abstol, dt) -> ev

Compute the characteristic multipliers of a continuous-time periodic matrix. 

For the given square periodic matrix `A(t)` of period `T`, 
the characteristic multipliers `ev` are the eigenvalues of 
the monodromy matrix `Ψ = Φ(T,0)`, where `Φ(t,τ)` is the state transition matrix satisfying the homogeneous linear ODE 

    dΦ(t,τ)/dt = A(t)Φ(t,τ),  Φ(τ,τ) = I. 

If `lifting = false`, `Ψ` is computed as a product of `K` state transition matrices 
`Ψ = Ψ_K*...*Ψ_1` (see [`monodromy`](@ref) with the associated keyword arguments). 
The eigenvalues are computed using the periodic Schur decomposition method of [1].

If `lifting = true`, `Ψ` is (implicitly) expressed as `Ψ = inv(N)*M`, where `M-λN` is a regular
pencil with `N` invertible and  
the eigenvalues of `M-λN` are the same as those of the matrix product
`Ψ := Ψ_K*...*Ψ_1`. 
An efficient version of the structure exploiting fast reduction method of [2] is employed, 
which embeds the determination of transition matrices into the reduction algorithm. 
This option may occasionally lead to inaccurate results for large values of `K`. 
`A` has the type [`FourierFunctionMatrix`](@ref).

_References_

[1] A. Bojanczyk, G. Golub, and P. Van Dooren, 
    The periodic Schur decomposition. Algorithms and applications, Proc. SPIE 1996.

[2] A. Varga & P. Van Dooren. Computing the zeros of periodic descriptor systems.
    Systems and Control Letters, 50:371-381, 2003.

"""
function pseig(at::PM, K::Int = 1; lifting::Bool = false, solver = "non-stiff", reltol = 1e-3, abstol = 1e-7, dt = at.period/100/at.nperiod) where 
   {T, PM <: FourierFunctionMatrix{:c,T}}
   n = size(at,1)
   n == size(at,2) || error("the periodic matrix must be square")
   nperiod = at.nperiod
   t = 0  
   Ts = at.period/K/nperiod
   if lifting 
      if K == 1
         ev = eigvals(tvstm(at, at.period, 0; solver, reltol, abstol, dt)) 
      else   
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
      end
      sorteigvals!(ev)
   else
      M = monodromy(at, K; solver, reltol, abstol, dt) 
      ev = K == 1 ? eigvals(view(M.M,:,:,1)) : pschur(M.M; withZ = false)[3]
      isreal(ev) && (ev = real(ev))
   end
   return nperiod == 1 ? ev : ev.^nperiod
end
function psceig(at::PM, K::Int = 1; kwargs...) where {T, PM <:FourierFunctionMatrix{:c,T}} 
   if isconstant(at)
      ce = eigvals(at(0))
   else
      ce = log.(complex(pseig(at, K; kwargs...)))/at.period
   end
   return isreal(ce) ? real(ce) : ce
end

"""
    psceigfr(A::FourierFunctionMatrix[, N]; P, atol) -> ce

Compute the characteristic exponents of a continuous-time periodic matrix in _Fourier Function Matrix_ representation. 

For a given square continuous-time periodic function matrix `A(t)` of period `T` 
in a  _Fourier Function Matrix_ representation, 
the characteristic exponents `ce` are computed as the eigenvalues of the state operator `A(t)-D*I` lying in the 
fundamental strip `-ω/2 <  Im(λ) ≤ ω/2`, where `ω = 2π/T`. A finite dimensional truncated matrix of order `n*(2*N*P+1)` 
is built to approximate `A(t)-D*I`, where `n` is the order of `A(t)`,  `N` is the number of selected harmonic components
in the Fourier representation and `P` is the period multiplication number (default: `P = 1`).
The default value used for `N` is `N = max(10,p-1)`, where `p` the number of harmonics terms of `A(t)` (see [`FourierFunctionMatrix`](@ref)). 

The keyword argument `atol` (default: `atol = 1.e-10`) is a tolerance on the magnitude of the trailing components of the 
associated eigenvectors used to validate their asymptotic (exponential) decay. Only eigenvalues satisfying this check are returned in `ce`. 
"""
function psceigfr(Afun::FourierFunctionMatrix{:c,T}, N::Int = max(10,maximum(ncoefficients.(Matrix(Afun.M)))); P::Int = 1, atol::Real = 1.e-10) where T
   n = size(Afun,1)
   n == size(Afun,2) || error("the periodic matrix must be square") 
   (N == 0 || isconstant(Afun)) && (return eigvals(getindex.(coefficients.(Matrix(Afun.M)),1)))
   Af = P == 1 ? Afun :  FourierFunctionMatrix(Fun(t -> Afun.M(t),Fourier(0..P*Afun.period)))
   D = Derivative(domain(Af.M))

   Aop = Af.M - DiagDerOp(D,n)
   NA = n*(2*N*P+1)
   RW = Aop[1:NA,1:NA]
   W = Matrix(RW)
   ev, V = eigen(W)

   ind = sortperm(imag(ev),by=abs) 
   ωhp2 = pi/Af.period/Af.nperiod
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
function DiagDerOp(D::Union{ApproxFunBase.DerivativeWrapper,ApproxFunBase.ConstantTimesOperator}, n::Int) 
   Z = tuple(D,ntuple(n->0I,n-1)...)
   for i = 2:n
       Z1 = tuple(ntuple(n->0I,i-1)...,D,ntuple(n->0I,n-i)...)
       Z = tuple(Z...,Z1...)
   end
   return hvcat(n,Z...)
end

# conversions
"""
     ffm2hr(Afun::FourierFunctionMatrix; atol = 0, rtol = √ϵ, squeeze = true) -> Ahr::HarmonicArray

Compute the harmonic (Fourier) representation of a Fourier periodic matrix object. 

The Fourier function matrix object `Afun` of period `T` is built using
the Fourier series representation of a periodic matrix `Afun(t)` of subperiod `T′ = T/k`, 
where each entry of `Afun(t)` has the form

             p
      a_0 +  ∑ ( ac_i*cos(i*t*2*π/T′)+as_i*sin(i*2*π*t/T′) ) ,
            i=1 

where `k ≥ 1` is the number of subperiods (default: `k = 1`).   

The harmonic array object `Ahr` of period `T` is built using
the harmonic representation of a periodic matrix `Ahr(t)` of subperiod `T′′ = T/k′` in the form

                     p′
     Ahr(t) = A_0 +  ∑ ( Ac_i*cos(i*t*2*π/T′′)+As_i*sin(i*2*π*t/T′′) ) ,
                    i=1 

where `p′` is the maximum index `j`, such that `Ac_j` and/or `As_j` are nonzero.
The tolerance used to assess nonzero elements is `tol = max(atol, rtol*maxnorm)`, where 
`maxnorm` is the maximum absolute value of the coefficients `ac_i` and `as_i` in `Afun(t)`. The default values of tolerances
are `atol = 0` and `rtol = √ϵ`, where `ϵ` is the working machine precision.
The resulting harmonic approximation `Ahr(t)` is returned in the harmonic array object `Ahr` 
(see [`HarmonicArray`](@ref)). 
"""
function ffm2hr(A::FourierFunctionMatrix{:c,T}; atol::Real = 0, rtol::Real = 0, squeeze::Bool = true) where  {T}
   lens = length.(coefficients.(Matrix(A.M)))
   #n = max(div(maximum(lens)-1,2),0)
   n = max(div(maximum(lens)-1,2),0)+1
   n1, n2 = size(A)

   ncur = n
   AHR = zeros(complex(T), n1, n2, n+1)
   tol = iszero(atol) ? (iszero(rtol) ? 10*n*eps(maximum(norm.(coefficients.(Matrix(A.M)),Inf))) : rtol*maximum(norm.(coefficients.(Matrix(A.M)),Inf)) ) : atol
   for i = 1:n1
       for j = 1:n2
           tt = coefficients(A.M[i,j])
           tt[abs.(tt) .<= tol] .= zero(T)
           lentt = lens[i,j]
           if lentt > 0
              AHR[i,j,1] = tt[1] 
              k = 1
              for k1 = 2:2:lentt-1
                  k += 1
                  AHR[i,j,k] = tt[k1+1] + im*tt[k1]
              end
              isodd(lentt) || (AHR[i,j,k+1] = im*tt[end])
           end
       end
   end
   nperiod = A.nperiod
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
      return HarmonicArray(AHR[:,:,[1;nperiod+1:nperiod:ncur]],A.period;nperiod)
   else
      return HarmonicArray(AHR[:,:,1:max(1,ncur+1)],A.period;nperiod)
   end       

end

"""
     tvmeval(A::FourierFunctionMatrix, t) -> Aval::Vector{Matrix}

Evaluate the time response of a periodic function matrix.

For the periodic matrix `A(t)`, in a Fourier Function Matrix representation, and the vector of time values `t`, 
`Aval[i] = A(t[i])` is evaluated for each time value `t[i]`. 
"""
function tvmeval(A::FourierFunctionMatrix, t::Union{Real,Vector{<:Real}} )
   te = isa(t,Real) ? [mod(t,A.period)] : mod.(t,A.period)
   return (A.M).(te)
end
function tpmeval(A::FourierFunctionMatrix, t::Real )
    (A.M)(t)
end
(F::FourierFunctionMatrix)(t) = (F.M)(t)
   
"""
    pmaverage(A::FourierFunctionMatrix) -> Am 

Compute for the continuous-time periodic matrix `A(t)` 
the corresponding time averaged matrix `Am` over one period.  
"""
function pmaverage(A::FourierFunctionMatrix)
   typeof(size(A.M)) == Tuple{} ? (return coefficients(A.M)[1]) : (return get.(coefficients.(Matrix(A.M)),1,0.0))
end
pmcopy(A::FourierFunctionMatrix) = FourierFunctionMatrix(A.M)

