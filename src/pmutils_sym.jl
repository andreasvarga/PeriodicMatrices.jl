function tvstm(A::PM, tf::Real, t0::Real = 0; solver = "", reltol = 1e-3, abstol = 1e-7, dt = (tf-t0)/10) where 
         {PM <: PeriodicSymbolicMatrix}
   n = size(A,1)
   n == size(A,2) || error("the function matrix must be square")

   isconstant(A) && ( return exp(tpmeval(A,t0)*(tf-t0)) )
   
   T1 = promote_type(typeof(t0), typeof(tf))
   T = Float64

   # using OrdinaryDiffEq
   u0 = Matrix{Float64}(I,n,n)
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
function monodromy(A::PM, K::Int = 1; solver = "non-stiff", reltol = 1e-3, abstol = 1e-7, dt = A.period/max(K,100)) where
         {PM <: PeriodicSymbolicMatrix} 
   n = size(A,1)
   n == size(A,2) || error("the periodic matrix must be square")
   nperiod = A.nperiod
   Ts = A.period/K/nperiod
   T = Float64

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
     pseigsm(PeriodicSymbolicMatrix[, K = 1]; lifting = false, solver, reltol, abstol, dt) -> ev

Compute the characteristic multipliers of a continuous-time periodic symbolic matrix. 

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
`A` has the type [`PeriodicSymbolicMatrix`](@ref).

_References_

[1] A. Bojanczyk, G. Golub, and P. Van Dooren, 
    The periodic Schur decomposition. Algorithms and applications, Proc. SPIE 1996.

[2] A. Varga & P. Van Dooren. Computing the zeros of periodic descriptor systems.
    Systems and Control Letters, 50:371-381, 2003.

"""
function pseigsm(at::PM, K::Int = 1; lifting::Bool = false, solver = "non-stiff", reltol = 1e-3, abstol = 1e-7, dt = at.period/100/at.nperiod) where 
   {PM <: PeriodicSymbolicMatrix}
   n = size(at,1)
   n == size(at,2) || error("the periodic matrix must be square")
   K > 0 || ArgumentError(throw("K must be positive: got K = $K"))
   nperiod = at.nperiod
   t = 0  
   Ts = at.period/K/nperiod
   T = Float64
   if K == 1
      ev = eigvals(tvstm(at, at.period, 0; solver, reltol, abstol, dt)) 
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
function pseig(at::PeriodicSymbolicMatrix, K::Int = 1; kwargs...) 
   pseig(convert(PeriodicFunctionMatrix,at),K; kwargs...)
end
"""
     psceigsm(A::PeriodicSymbolicMatrix[, K = 1]; lifting = false, solver, reltol, abstol, dt) -> ce

Compute the characteristic exponents of a periodic symbolic matrix.

For a given square continuous-time periodic function matrix `A(t)` of period `T`, 
the characteristic exponents `ce` are computed as `log.(ev)/T`, 
where  `ev` are the characteristic
multipliers (i.e., the eigenvalues of the monodromy matrix of `A(t)`).  
For available options see [`pseig(::PeriodicFunctionMatrix)`](@ref). 
"""
function psceigsm(at::PeriodicSymbolicMatrix, K::Int = 1; kwargs...) 
   ce = log.(complex(pseigsm(at, K; kwargs...)))/at.period
   return isreal(ce) ? real(ce) : ce
end
function psceig(at::PeriodicSymbolicMatrix, K::Int = 1; kwargs...) 
   ce = log.(complex(pseig(at, K; kwargs...)))/at.period
   return isreal(ce) ? real(ce) : ce
end

# conversions
"""
     psm2hr(A::PeriodicSymbolicMatrix; nsample, NyquistFreq) -> Ahr::HarmonicArray

Convert a periodic symbolic matrix into a harmonic array representation. 
If `A(t)` is a periodic symbolic matrix of period `T`, then the harmonic array representation
`Ahr` is determined by applying the fast Fourier transform to the sampled arrays `A(iΔ)`, `i = 0, ..., k`,
where `Δ = T/k` is the sampling period and `k` is the number of samples specified by the keyword argument
`nsample = k` (default: `k = 128`). If the Nyquist frequency `f` is specified via the keyword argument
`NyquistFreq = f`, then `k` is chosen `k = 2*f*T` to avoid signal aliasing.     
"""
function psm2hr(A::PeriodicSymbolicMatrix; nsample::Int = 128, NyquistFreq::Union{Real,Missing} = missing)   
   isconstant(A) && (return HarmonicArray(convert(PeriodicFunctionMatrix,A).f(0),A.period; nperiod = A.nperiod))
   nsample > 0 || ArgumentError("nsample must be positive, got $nsaple")
   ns = ismissing(NyquistFreq) ? nsample : Int(floor(2*abs(NyquistFreq)*A.period))+1
   Δ = A.period/ns
   ts = (0:ns-1)*Δ
   return ts2hr(PeriodicTimeSeriesMatrix(convert(PeriodicFunctionMatrix,A).f.(ts), A.period; nperiod = A.nperiod))
end
"""
     tvmeval(Asym::PeriodicSymbolicMatrix, t) -> A::Vector{Matrix}

Evaluate the time response of a periodic symbolic matrix.

For the periodic symbolic matrix `Asym` representing a time periodic matrix `A(t)`
and the vector of time values `t`, 
`A[i] = A(t[i])` is evaluated for each time value `t[i]`. 
"""
function tvmeval(A::PeriodicSymbolicMatrix, t::Union{Real,Vector{<:Real}} )
   te = isa(t,Real) ? [mod(t,A.period)] : mod.(t,A.period)
   return (convert(PeriodicFunctionMatrix,A).f).(te)
end
function tpmeval(A::PeriodicSymbolicMatrix, t::Real)
   te = isa(t,Real) ? [mod(t,A.period)] : mod.(t,A.period)
   return (convert(PeriodicFunctionMatrix,A).f).(te)[1]
end
(F::PeriodicSymbolicMatrix)(t) = tpmeval(F, t)
   
"""
    pmaverage(A::PeriodicSymbolicMatrix) -> Am 

Compute for the continuous-time periodic matrix `A(t)` 
the corresponding time averaged matrix `Am` over one period.  
"""
function pmaverage(A::PM; rtol = sqrt(eps())) where {PM <: PeriodicSymbolicMatrix} 
   tsub = A.period/A.nperiod
   tt, = quadgk(t -> tpmeval(A,t), 0., tsub; rtol)
   return tt/tsub
end
# function pmaverage(A::PeriodicSymbolicMatrix) 
#    return real(convert(HarmonicArray,A).values[:,:,1])
# end

"""
     hreval(Ahr::HarmonicArray, t; ntrunc, exact = true) -> A::Matrix

Evaluate the harmonic array `Ahr` representing a continuous-time 
time periodic matrix `A(t)` for a symbolic time value `t`. 
A symbolic evaluation of `A(t)` is performed (see also [`hr2psm`](@ref))
The keyword argument `ntrunc` specifies the number of harmonics to be used for the evaluation 
(default: maximum possible number). 
"""
function PeriodicMatrices.hreval(ahr::HarmonicArray{:c,T}, t::Num; ntrunc::Int = max(size(ahr.values,3)-1,0)) where {T}
      (ntrunc < 0 || ntrunc >= size(ahr.values,3)) && error("ntrunc out of allowed range")
   return hr2psm(ahr, 0:ntrunc)
end   

"""
     hr2psm(Ahr::HarmonicArray, nrange) -> A::Matrix{Num}

Convert a range of harmonic components `nrange` of the harmonic array `Ahr` to a symbolic matrix `A`. 
The default range is `nrange = 0:n`, where `n` is the order of the maximum harmonics.  
"""
function hr2psm(ahr::HarmonicArray, nrange::UnitRange = 0:size(ahr.values,3)-1)   
   Symbolics.@variables t   
   Period = ahr.period
   i1 = max(first(nrange),0)
   i2 = min(last(nrange), size(ahr.values,3)-1)
   
   ts = t*2*pi*ahr.nperiod/Period
   
   i1 == 0 ? a = Num.(real.(ahr.values[:,:,1])) : (i1 -=1; a = zeros(Num,size(ahr.values,1),size(ahr.values,2)))
   for i in i1+1:i2
       ta = view(ahr.values,:,:,i+1)
       a .+= real.(ta).*(cos(i*ts)) .+ imag.(ta) .* (sin(i*ts))
   end
   return a
end 

