module Test_pmutils

using ApproxFun
using Symbolics
using PeriodicMatrices
using Test
using LinearAlgebra
using MatrixPencils
using Bessels
#using BenchmarkTools

println("Test_pmutils")

@testset "test_pmutils" begin


# SwitchingPeriodicMatrix      
n = 2; pa = 3; px = 6; T = 10; 
Ad = 0.5*SwitchingPeriodicMatrix([rand(Float64,n,n) for i in 1:pa],[10,15,20],T);
@test pseig(Ad) ≈ pseig(Ad,fast=true)
@test pseig(convert(PeriodicArray,Ad).M,rev = false) ≈ pseig(convert(PeriodicArray,Ad).M,fast=true, rev = false)
ad = rand(2,2); A1d = PeriodicArray(ad,2)
@test psceig(A1d) ≈ eigvals(ad)

# Periodic Matrix with time-varying dimensions 
na = [5, 3, 3, 4, 1]; ma = [3, 3, 4, 1, 5]; pa = 5; px = 5;   
#na = 5*na; ma = 5*ma;
Ad = PeriodicMatrix([rand(Float64,ma[i],na[i]) for i in 1:pa],pa); 
@test  pseig(pm2pa(Ad)) ≈ pseig(Ad)

A = [rand(2,2), rand(2,2)]; B = [rand(2,2), rand(2,2)]
Ar1, Br1 = PeriodicMatrices.psreduc_fast(A,B)
Ar2, Br2 = PeriodicMatrices.psreduc_reg(A,B)
Ar3, Br3 = PeriodicMatrices.psreduc_reg(reshape(hcat(A...),2,2,2),reshape(hcat(B...),2,2,2))
ev = sort(eigvals(Ar1,Br1),by=abs)[1:2]
ev2 = eigvals(Ar2,Br2)
ev3 = eigvals(Ar3,Br3)
@test sort(real(ev)) ≈ sort(real(ev2)) ≈ sort(real(ev3))
@test sort(imag(ev)) ≈ sort(imag(ev2)) ≈ sort(imag(ev3))

Ad = reshape(vcat(A...),2,2,2); Bd = reshape(vcat(B...),2,2,2); 

# symbolic periodic 
@variables t
A = [cos(t) 1; 1 1-sin(t)];
B = [cos(t)+sin(t); 1-sin(t)];
C = [sin(t)+cos(2*t) 1];
As = PeriodicSymbolicMatrix(A,2*pi);
Bs = PeriodicSymbolicMatrix(B,2*pi);
Cs = PeriodicSymbolicMatrix(C,2*pi);

# functional expressions
tA(t::Real) = [cos(t) 1; 1 1-sin(t)];
tB(t::Real) = [cos(t)+sin(t); 1-sin(t)];
tC(t::Real) = [sin(t)+cos(2*t) 1];

At = PeriodicFunctionMatrix(tA,2pi)
Bt = PeriodicFunctionMatrix(tB,2pi)
Ct = PeriodicFunctionMatrix(tC,2pi)

s = Fourier(0..2π)
Af = FourierFunctionMatrix(Fun(tA,s), 2pi)


# store snapshots as 3d arrays
N = 200; 
tg = collect((0:N-1)*2*pi/N);
ts = (0:N-1)*2*pi/N;
# At = reshape(hcat(A.(t)...),2,2,N);  
# Bts = reshape(hcat(B.(t)...),2,1,N);  
# Ct = reshape(hcat(C.(t)...),1,2,N);

# time series expressions
Ats = PeriodicTimeSeriesMatrix(tA.(ts),2*pi);
Bts = PeriodicTimeSeriesMatrix(tB.(ts),2*pi);
Cts = PeriodicTimeSeriesMatrix(tC.(ts),2*pi);

# harmonic expressions
@time Ahr = ts2hr(Ats);
@test Ahr.values[:,:,1] ≈ [0. 1; 1 1] && Ahr.values[:,:,2] ≈ [1. 0; 0 -im] 
@time Bhr = ts2hr(Bts);
@test Bhr.values[:,:,1] ≈ [0.; 1] && Bhr.values[:,:,2] ≈ [1.0+im; -im]
@time Chr = ts2hr(Cts);
@test Chr.values[:,:,1] ≈ [0. 1] && Chr.values[:,:,2] ≈ [im 0] && Chr.values[:,:,3] ≈ [1 0]

#@time Affm = ts2ffm(Ats);

nperiod = 24
time1 = (0:N-1)*2*pi*nperiod/N;
Ats1 = PeriodicTimeSeriesMatrix(tA.(time1),2*pi*nperiod);
Ahr1 = ts2hr(Ats1);
@test convert(PeriodicFunctionMatrix,Ahr1).f(1) ≈ convert(PeriodicFunctionMatrix,Ahr).f(1) 

@test iszero(hr2psm(Ahr,1:1) + hr2psm(Ahr,0:0) - hr2psm(Ahr))
@test iszero(hr2psm(Ahr1,1:1) + hr2psm(Ahr1,0:0) - hr2psm(Ahr1))

# harmonic vs. symbolic
@test norm(substitute.(convert(PeriodicSymbolicMatrix,Ahr).F - A, (Dict(t => rand()),))) < 1e-15
@test norm(substitute.(convert(PeriodicSymbolicMatrix,Bhr).F - B, (Dict(t => rand()),))) < 1e-15
@test norm(substitute.(convert(PeriodicSymbolicMatrix,Chr).F - C, (Dict(t => rand()),))) < 1e-15

# harmonic vs. time series
@test all(norm.(tvmeval(Ahr,tg).-Ats.values) .< 1.e-7)
@test all(norm.(tvmeval(Bhr,tg).-Bts.values) .< 1.e-7)
@test all(norm.(tvmeval(Chr,tg).-Cts.values) .< 1.e-7)

# check time values on the grid
for method in ("constant", "linear", "quadratic", "cubic")
      @test all(norm.(tvmeval(Ats,tg; method).-Ats.values) .< 1.e-7)
      @test all(norm.(tvmeval(Bts,tg; method).-Bts.values) .< 1.e-7)
      @test all(norm.(tvmeval(Cts,tg; method).-Cts.values) .< 1.e-7)
end   

# check interpolated values: time series vs. harmonic
tt = rand(10)*2pi;
for method in ("linear", "quadratic", "cubic")
    @test all(norm.(tvmeval(Ats,tt; method).-tvmeval(Ahr,tt; exact = true)) .< 1.e-3)
    @test all(norm.(tvmeval(Bts,tt; method).-tvmeval(Bhr,tt; exact = false)) .< 1.e-3)
    @test all(norm.(tvmeval(Cts,tt; method).-tvmeval(Chr,tt; exact = true)) .< 1.e-3)
end

# check conversion to function form
Amat = convert(PeriodicFunctionMatrix,Ahr); 
@test all(norm.(tvmeval(Ats,tt; method = "linear").-tvmeval(Ahr,tt)) .< 1.e-3)
@test iszero(convert(PeriodicSymbolicMatrix,Amat).F-As.F)

Amat =  convert(PeriodicFunctionMatrix,Ats);
@test all(norm.(tvmeval(Ats,tt; method = "linear").-Amat.f.(tt)) .< 1.e-3)
#@test iszero(convert(PeriodicSymbolicMatrix,Amat).F-As.F)

Amat = PeriodicFunctionMatrix(tA,2pi);
@test all(norm.(tvmeval(Ats,tt; method = "linear").-tvmeval(Amat,tt)) .< 1.e-3)
@test iszero(convert(PeriodicSymbolicMatrix,Amat).F-As.F)

Amat =  convert(PeriodicFunctionMatrix,As);
@test all(norm.(tvmeval(Ats,tt; method = "linear").-tvmeval(As,tt)) .< 1.e-3)
@test iszero(convert(PeriodicSymbolicMatrix,Amat).F-As.F)
@test size(Amat) == size(As)


for method in ("constant", "linear", "quadratic", "cubic")
      Amat = ts2pfm(Ats; method);
      @test all(norm.(Ats.values.-Amat.f.(tg)) .< 1.e-10)
      Bmat = ts2pfm(Bts; method);
      @test all(norm.(Bts.values.-Bmat.f.(tg)) .< 1.e-10)
      Cmat = ts2pfm(Cts; method);
      @test all(norm.(Cts.values.-Cmat.f.(tg)) .< 1.e-10)
end
@test_throws ArgumentError ts2pfm(Ats,method = "noidea")

Ats1 = [Ats.values;[Ats.values[1]]]; N = length(Ats1); 
d = 2pi/(N-1); ts1 = (0:N-1)*d
At0=PeriodicMatrices.ts2fm(Ats1,2pi,method = "constant");
@test At0.(ts1) ≈ (PeriodicMatrices.ts2fm(Ats1,2pi,method = "linear")).(ts1);
@test At0.(ts1) ≈ (PeriodicMatrices.ts2fm(Ats1,2pi,method = "quadratic")).(ts1);
@test At0.(ts1) ≈ (PeriodicMatrices.ts2fm(Ats1,2pi,method = "cubic")).(ts1);
@test_throws ArgumentError PeriodicMatrices.ts2fm(Ats1,2pi,method = "noidea")

# test pseig
@test sort(pseig(At,200)) ≈ sort(pseig(As,200)) ≈ sort(pseig(Ahr,200)) ≈ sort(pseig(Af,200))
@test sort(pseig(At,1)) ≈ sort(pseig(As,1)) ≈ sort(pseig(Ahr,1)) ≈ sort(pseig(Af,1))
@test sort(pseig(At,200,lifting = true)) ≈ sort(pseig(As,200,lifting = true)) ≈ sort(pseig(Ahr,200,lifting = true)) ≈ sort(pseig(Af,200,lifting = true))
@test norm(pseig(convert(PeriodicFunctionMatrix,Ats;method="constant"),length(Ats); abstol = 1.e-14,reltol=1.e-14)-pseig(Ats)) < 1.e-7
@test norm(sort(psceig(At,200,abstol=1.e-14,reltol = 1.e-14))-sort(psceig(Ats)),Inf) < 1.e-5

@test sort(psceig(At,200)) ≈ sort(psceig(As,200)) ≈ sort(psceig(Ahr,200)) ≈ sort(psceig(Af,200))

aa = rand(2,2); ah = HarmonicArray(aa,2); @test psceig(ah) ≈ eigvals(aa)
@test norm(sort(psceig(Ahr,200)) - sort(psceighr(Ahr,20),by=real)) < 1.e-7

@test ts2hr(Ats;n=0)(0) ≈ pmaverage(Ats)

@test hr2bt(Ahr,4)[1:10,1:10] == hr2bt(Ahr,2)
@test hr2btupd(Ahr,4) == hr2btupd(Ahr,4)


# example of Colaneri
at(t) = [0 1; -10*cos(t) -24-10*sin(t)];
Afun=PeriodicFunctionMatrix(at,2pi); 
ev = pseig(Afun; solver = "non-stiff", reltol = 1.e-10, abstol = 1.e-10)
# @time Φ = tvstm(Afun.f, Afun.period; reltol = 1.e-10)
cvals = log.(complex(ev))/2/pi 
println("cvals = $cvals  -- No one digit accuracy!!!")
@test maximum(abs.(ev)) ≈ 1

# using ApproxFun
s = Fourier(0..2π)
Af = FourierFunctionMatrix(Fun(t -> [0 1; -10*cos(t) -24-10*sin(t)],s), 2pi)
Atfun = convert(PeriodicFunctionMatrix,Af)
Ahrfun = convert(PeriodicFunctionMatrix,pfm2hr(Afun))


@time cvals0 = psceig(Af, 500; solver = "non-stiff", reltol = 1.e-10, abstol = 1.e-10)
@time cvals = psceig(Afun, 500; solver = "non-stiff", reltol = 1.e-10, abstol = 1.e-10)
@time cvals1 = psceig(Atfun, 500; solver = "non-stiff", reltol = 1.e-10, abstol = 1.e-10)
@time cvals2 = psceig(Ahrfun, 500; solver = "non-stiff", reltol = 1.e-10, abstol = 1.e-10)
@time cvals4 = psceigfr(Af,40)
@test sort(cvals0) ≈ sort(cvals) ≈ sort(cvals1) ≈ sort(cvals2) && norm(sort(cvals4)-sort(cvals0)) < 1.e-6

solver = "non-stiff"
#for solver in ("non-stiff", "stiff", "linear", "symplectic", "noidea")
for solver in ("non-stiff", "stiff", "linear", "noidea")
      println("solver = $solver")
    @time cvals = psceig(Af, 500; solver, reltol = 1.e-10, abstol = 1.e-10)
    @test isapprox(cvals, [0; -24], atol = 1.e-7)
end

Af = FourierFunctionMatrix(Fun(t -> [0 1; 0. -24],s), 2pi)
cvals = psceig(Af)
@test isapprox(sort(cvals), [-24.; 0.0], atol = 1.e-7)
evals = pseig(Af)
@test isapprox(sort(evals), sort(eigvals(exp(Af(0)*2*pi))), atol = 1.e-60)

# Tt = Fun(t -> [12+5*sin(t) 1/2; 1 0],s)
# Tinvt=inv(Tt)
# Atilde=Tt*At.M*Tinvt+Tt'*Tinvt
# Aref = Fun(t -> [0 0; 2 -24-10*sin(t)],s)
# @test norm(Aref-Atilde) < 1.e-10

# example Floquet analysis from ApproxFun.jl
a=0.15
at1(t) = -[0 -1 0 0; (2+a*cos(2t)) 0 -1 0; 0 0 0 -1; -1 0 (2+a*cos(2t)) 0]
Afun1=PeriodicFunctionMatrix(at1,pi);
ev1 = pseig(Afun1; solver = "non-stiff", reltol = 1.e-10, abstol = 1.e-10)
cvals1 = log.(complex(ev1))/pi

Af = FourierFunctionMatrix(Fun(at1,Fourier(0..π)), pi)
ev2 = pseig(Af,200,lifting = false)
ev3 = pseig(Af,200,lifting = true)
@test sort(real(ev1)) ≈ sort(real(ev2)) ≈ sort(real(ev3))
@test sort(imag(ev1)) ≈ sort(imag(ev2)) ≈ sort(imag(ev3))

a=0.15
at2(t) = -[0 -1 0 0; (2+a*cos(2t)) 0 -1 0; 0 0 0 -1; -1 0 (2+a*cos(2t)) 0]
Afun2=PeriodicFunctionMatrix(at2,2*pi;nperiod=2);
ev2 = pseig(Afun2; solver = "non-stiff", reltol = 1.e-10, abstol = 1.e-10)
cvals2 = log.(complex(ev2))/(2pi)
@test ev1.^2 ≈ ev2 && real(cvals1) ≈ real(cvals2) && Afun1 !== Afun2


# full accuracy characteristic exponents
# solver = "symplectic"
# for solver in ("non-stiff", "stiff", "linear", "symplectic", "noidea")
#     @time M = monodromy(Afun, 500; solver, reltol = 1.e-10, abstol = 1.e-10);
#     cvals = log.(complex(pseig(M)))/2/pi
#     #println("solver = $solver cvals = $cvals")
#     @test isapprox(cvals, [0; -24], atol = 1.e-7)
# end

solver = "non-stiff"
#for solver in ("non-stiff", "stiff", "linear", "symplectic", "noidea")
for solver in ("non-stiff", "stiff", "linear", "noidea")
    println("solver = $solver")
    @time cvals = psceig(Afun, 500; solver, reltol = 1.e-10, abstol = 1.e-10)
    @test isapprox(cvals, [0; -24], atol = 1.e-7)
end

# Vinograd example: unstable periodic system with all A(t) stable
#at(t) = [3.5 6;-6 -5.5]+[-4.5 6; 6 4.5]*cos(12*t)+[6 4.5;4.5 -6]*sin(12*t); T = pi/6; # does not work
function at(t::Real)
    [-1-9*(cos(6*t))^2+12*sin(6*t)*cos(6*t) 12*(cos(6*t))^2+9*sin(6*t)*cos(6*t);
    -12*(sin(6*t))^2+9*sin(6*t)*cos(6*t) -1-9*(sin(6*t))^2-12*sin(6*t)*cos(6*t)]
end
@test eigvals(at(rand())) ≈ [-10,-1]
T = pi/3;
Afun=PeriodicFunctionMatrix(at,T);
#for solver in ("non-stiff", "stiff", "linear", "symplectic", "noidea")
solver = "linear"
for solver in ("non-stiff", "stiff", "linear", "noidea")
      println("solver = $solver")
      @time cvals = psceig(Afun, 500; solver, reltol = 1.e-10)
      @test cvals ≈ [2; -13]
end  

@variables t
A11 =  [-1-9*(cos(6*t))^2+12*sin(6*t)*cos(6*t) 12*(cos(6*t))^2+9*sin(6*t)*cos(6*t);
      -12*(sin(6*t))^2+9*sin(6*t)*cos(6*t) -1-9*(sin(6*t))^2-12*sin(6*t)*cos(6*t)]

As = PeriodicSymbolicMatrix(A11,pi/3)      

@test eigvals(As(rand())) ≈ [-10,-1]
@test pmaverage(As) ≈ pmaverage(Afun)

  #for solver in ("non-stiff", "stiff", "linear", "symplectic", "noidea")
solver = "linear"
for solver in ("non-stiff", "stiff", "linear", "noidea")
        println("solver = $solver")
        @time cvals = psceig(As, 500; solver, reltol = 1.e-10)
        @test sort(cvals) ≈ sort([2; -13])
        @time cvals = psceigsm(As, 500; solver, reltol = 1.e-10)
        @test sort(cvals) ≈ sort([2; -13])
        @time cvals = psceigsm(As, 500; lifting = true, solver, reltol = 1.e-10)
        @test sort(cvals) ≈ sort([2; -13])
end 
@time cvals = psceigsm(As, 1; solver, reltol = 1.e-10)
@test ≈(sort(cvals),sort([2; -13]),atol=1.e-6)


  

# Examples from Richards's Book
a = 1; q = .1; τ = pi/3;

ψ1(t) = t <= τ ? t/τ : 0. # Hill equation 
ψ2(t) = cos(2t)  # Mathieu equation

# periodic function 
A1(t) = [0. 1.; -a+2*q*ψ1(t) 0]
A2(t) = [0. 1.; -a+2*q*ψ2(t) 0 ]

G1f = PeriodicFunctionMatrix(A1,pi)
cm1f = pseig(G1f,1000; reltol = 1.e-14)

G2f = PeriodicFunctionMatrix(A2,pi)
cm2f = pseig(G2f,1000; reltol = 1.e-14)


# Floquet analysis for Meissner equation of Example Fig 3.1
a = 1; q = .1; τ = pi/3;

# periodic function based approach
ψ3(t) = t <= τ ? 1. : -1. # Meissner equation Example Fig 3.1
A3(t) = [0. 1.; -a+2*q*ψ3(t) 0]
G3f = PeriodicFunctionMatrix(A3,pi)
cm3f = pseig(G3f,1000; reltol = 1.e-14) # characteristic multipliers

# analytic approach form Richards's book page 29
c = sqrt(a-2q); d = sqrt(a+2q)
Φ = [cos(d*(pi-τ)) 1/d*sin(d*(pi-τ)); -d*sin(d*(pi-τ)) cos(d*(pi-τ)) ]*[cos(c*τ) 1/c*sin(c*τ); -c*sin(c*τ) cos(c*τ)]
cm3 = eigvals(Φ)  # exact characteristic multipliers
ce3 = log.(complex(cm3))/pi

# periodic switching matrix based approach
G3sw = PeriodicSwitchingMatrix([[0 1; -a+2q 0],[0 1; -a-2q 0]],[0,τ],pi)
cm3sw = pseig(G3sw) # characteristic multipliers 
ce3sw = psceig(G3sw) # characteristic exponents 
@test norm(sort(cm3f)-sort(cm3),Inf) < 1.e-6 && norm(sort(cm3sw)-sort(cm3),Inf) < 1.e-14 && norm(sort(ce3sw,by=real)-sort(ce3,by=real),Inf) < 1.e-14


# Floquet analysis for Hill equation with a negative slope sawtooth waveform coefficient of Example Fig 3.3
a = 1; q = .1; 

# periodic function based approach
ψ4(t) = -2*t/pi+1 # Hill equation Example Fig 3.3
A4(t) = [0. 1.; -a+2*q*ψ4(t) 0]
G4f = PeriodicFunctionMatrix(A4,pi)
cm4f = pseig(G4f,1000; reltol = 1.e-14) # characteristic multipliers


# analytic approach form Richards's book page 35
c = sqrt(a-2q); d = sqrt(a+2q); k = pi/6/q
σd = k*d^3; σc = k*c^3
Φ4 = (pi^2/(6*sqrt(3)*q))*[d*(a-2q)*(besselj(1/3,σd)*besselj(2/3,σc)+besselj(-1/3,σd)*besselj(-2/3,σc)) c*d*(besselj(1/3,σd)*besselj(-1/3,σc)-besselj(-1/3,σd)*besselj(1/3,σc));
                        (a+2q)*(a-2q)*(besselj(-2/3,σd)*besselj(2/3,σc)-besselj(2/3,σd)*besselj(-2/3,σc)) c*(a+2q)*(besselj(-2/3,σd)*besselj(-1/3,σc)+besselj(2/3,σd)*besselj(1/3,σc))]
cm4 = eigvals(Φ4)  # exact characteristic multipliers
@test norm(sort(cm4f)-sort(cm4),Inf) < 1.e-7

# Floquet analysis for Hill equation with a positive slope sawtooth waveform coefficient of Example Fig 3.4
a = 1; q = .1; 

# periodic function based approach
ψ5(t) = 2*t/pi-1 # Hill equation Example Fig 3.3
A5(t) = [0. 1.; -a+2*q*ψ5(t) 0]
G5f = PeriodicFunctionMatrix(A5,pi)
cm5f = pseig(G5f,1000; reltol = 1.e-14) # characteristic multipliers


# analytic approach form Richards's book page 36
using Bessels
c = sqrt(a-2q); d = sqrt(a+2q); k = pi/6/q
σd = k*d^3; σc = k*c^3
Φ5 = (-pi^2/(6*sqrt(3)*q))*[-d*(a-2q)*(besselj(1/3,σd)*besselj(2/3,σc)+besselj(-1/3,σd)*besselj(-2/3,σc)) c*d*(besselj(1/3,σd)*besselj(-1/3,σc)-besselj(-1/3,σd)*besselj(1/3,σc));
                        (a+2q)*(a-2q)*(besselj(-2/3,σd)*besselj(2/3,σc)-besselj(2/3,σd)*besselj(-2/3,σc)) -c*(a+2q)*(besselj(-2/3,σd)*besselj(-1/3,σc)+besselj(2/3,σd)*besselj(1/3,σc))]
cm5 = eigvals(Φ5)  # exact characteristic multipliers
@test norm(sort(cm5f)-sort(cm5),Inf) < 1.e-7


end
end

