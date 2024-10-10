module Test_pmops

using ApproxFun
using Symbolics
using PeriodicMatrices
using Test
using LinearAlgebra

println("Test_pmops")

@testset "pmops" begin


# generate periodic function matrices with same period and nperiod = 1
A(t) = [0  1; -10*cos(t)-1 -24-19*sin(t)]
X(t) = [1+cos(t) 0; 0 1+sin(t)]  # desired solution
Xder(t) = [-sin(t) 0; 0 cos(t)]  # derivative of the desired solution
C(t) = [ -sin(t)  -1-sin(t)-(-1-10cos(t))*(1+cos(t));
-1-sin(t)-(-1-10cos(t))*(1+cos(t))   cos(t)- 2(-24 - 19sin(t))*(1 + sin(t)) ]  # corresponding C
Cd(t) = [ sin(t)  -1-cos(t)-(-1-10cos(t))*(1+sin(t));
-1-cos(t)-(-1-10cos(t))*(1+sin(t))    -cos(t)-2(-24-19sin(t))*(1 + sin(t)) ] # corresponding Cd
A1(t) = [0  1; -0.5*cos(t)-1 -24-19*sin(t)]  # invertible matrix
A1inv(t) = [(24 + 19sin(t)) / (-1 - 0.5cos(t))  1 / (-1 - 0.5cos(t)); 1.0  0.0]  # invertible matrix

a(t) = cos(t)
@test PeriodicFunctionMatrix{:c,BigFloat}(a,2pi)(1) ≈ PeriodicFunctionMatrix(a,2pi)(1)
b(t) = [cos(t)]
@test PeriodicFunctionMatrix{:c,BigFloat}(b,2pi)(1) ≈ PeriodicFunctionMatrix(b,2pi)(1)

# PeriodicFunctionMatrix
@show "periodic function matrix operations"
At = PeriodicFunctionMatrix(A,2*pi)
Ct = PeriodicFunctionMatrix(C,2*pi)
Cdt = PeriodicFunctionMatrix(Cd,2*pi)
Xt = PeriodicFunctionMatrix(X,2*pi)
Xdert = PeriodicFunctionMatrix(Xder,2*pi)

@test eltype(At) == Float64
@test convert(PeriodicFunctionMatrix{:c,BigFloat},At)(1) ≈ At(1)
@test convert(HarmonicArray{:c,Float64},At) == convert(HarmonicArray,At)
@test convert(PeriodicSwitchingMatrix,At)(0) ≈ At(0)

ap = rand(2,2); Ap = PeriodicFunctionMatrix(ap,pi)
@test promote_period(ap,At) == promote_period(At) == At.period
@test promote_period(ap,Ap,At) == At.period
@test promote_period2(ap,At) == promote_period2(At) == (At.period, At.nperiod)
@test promote_period2(ap,Ap,At) == (At.period,At.nperiod)



@test set_period(set_period(At,4pi),2pi) == At
@test PeriodicFunctionMatrix(PeriodicFunctionMatrix(At,4pi),2pi) == At

@test iszero(PeriodicFunctionMatrix{:c,Float64}(zeros(Int,2,2),2))
@test eltype(At) == Float64
@test At[:,:] == At

@test At*Xt+Xt*At'+Ct ≈  pmderiv(Xt) ≈ Xdert
@test At'*Xt+Xt*At+Cdt ≈ -pmderiv(Xt)
@test norm(At*Xt+Xt*At'+Ct-pmderiv(Xt)) < 1.e-7
@test norm(At'*Xt+Xt*At+Cdt+pmderiv(Xt)) < 1.e-7



At = PeriodicFunctionMatrix(A,4*pi,nperiod=2)
Ct = PeriodicFunctionMatrix(C,2*pi)
Cdt = PeriodicFunctionMatrix(Cd,2*pi)
Xt = PeriodicFunctionMatrix(X,8*pi,nperiod=4)
Xdert = PeriodicFunctionMatrix(Xder,8*pi,nperiod=4)
@test At*Xt+Xt*At'+Ct ≈  pmderiv(Xt) ≈ Xdert
@test At'*Xt+Xt*At+Cdt ≈ -pmderiv(Xt)
@test norm(At*Xt+Xt*At'+Ct-pmderiv(Xt)) < 1.e-7
@test norm(At'*Xt+Xt*At+Cdt+pmderiv(Xt)) < 1.e-7

t = rand(); 
@test [At Ct](t) ≈ [At(t) Ct(t)]
@test [At; Ct](t) ≈ [At(t); Ct(t)]
@test blockdiag(At,Ct)(t) ≈ bldiag(At(t),Ct(t))


D = rand(2,2)
@test At+I == I+At && At*5 == 5*At && At*D ≈ -At*(-D) && iszero(At-At) && !iszero(At)
@test inv(At)*At ≈ I ≈ At*inv(At) && At+I == I+At
@test norm(At-At,1) == norm(At-At,2) == norm(At-At,Inf) == 0
@test iszero(opnorm(At-At,1)) && iszero(opnorm(At-At,2)) && iszero(opnorm(At-At,Inf)) && iszero(opnorm(At-At))
@test trace(At-At) == 0 && iszero(tr(At-At))

@test PeriodicFunctionMatrix(D,2*pi) !== PeriodicFunctionMatrix(D,4*pi) && 
      !(PeriodicFunctionMatrix(D,2*pi) ≈ PeriodicFunctionMatrix(D,4*pi))

@test tpmeval(At,1)[1:2,1:1] == tpmeval(At[1:2,1],1) && lastindex(At,1) == 2 && lastindex(At,2) == 2

# HarmonicArray
@show "Harmonic array operations"
A0 = rand(2,2); Acos = [rand(2,2)]; Asin = [rand(2,2),rand(2,2)]
HarmonicArray(A0,Acos,Asin,pi)
@test iszero(imag(HarmonicArray(A0,Acos,pi).values))
@test iszero(real(HarmonicArray(zeros(2,2),nothing,Asin,pi).values))
@test HarmonicArray(A0,pi) == HarmonicArray(A0,nothing,nothing,pi) == HarmonicArray{:c,Float64}(A0,pi)

At = PeriodicFunctionMatrix(A,2*pi); 
Ah = convert(HarmonicArray,At);
@test convert(PeriodicFunctionMatrix{:c,BigFloat},Ah)(1) ≈ Ah(1)
@test convert(PeriodicTimeSeriesMatrix,Ah) ≈ convert(PeriodicTimeSeriesMatrix,At)


Ah1 = convert(HarmonicArray,PeriodicFunctionMatrix(A1,2*pi));
Ch = convert(HarmonicArray,PeriodicFunctionMatrix(C,2*pi));
Cdh = convert(HarmonicArray,PeriodicFunctionMatrix(Cd,2*pi));
Xh = convert(HarmonicArray,PeriodicFunctionMatrix(X,2*pi));
Xderh = convert(HarmonicArray,PeriodicFunctionMatrix(Xder,2*pi));

@test set_period(set_period(Ah,4pi),2pi) == Ah
@test PeriodicFunctionMatrix(PeriodicFunctionMatrix(At,4pi),2pi) == At

@test iszero(PeriodicFunctionMatrix{:c,Float64}(zeros(Int,2,2),2))
@test eltype(At) == Float64
@test At[:,:] == At

@test issymmetric(Ch) && issymmetric(Cdh) && issymmetric(Xh) && issymmetric(Xderh)
@test Ah*Xh+Xh*Ah'+Ch ≈  pmderiv(Xh) ≈ Xderh
@test Ah'*Xh+Xh*Ah+Cdh ≈ -pmderiv(Xh) 
@test norm(Ah*Xh+Xh*Ah'+Ch-pmderiv(Xh),Inf) < 1.e-7 && norm(pmderiv(Xh)- Xderh) < 1.e-7
@test norm(Ah'*Xh+Xh*Ah+Cdh+pmderiv(Xh),1) < 1.e-7

D = rand(2,2)
@test Ah+I == I+Ah && Ah*5 == 5*Ah && Ah*D ≈ -Ah*(-D) && iszero(Ah-Ah) && !iszero(Ah)
@test HarmonicArray(D,2*pi) !== HarmonicArray(D,4*pi) && 
      !(HarmonicArray(D,2*pi) ≈ HarmonicArray(D,4*pi))
Ah1i = inv(Ah1)
@test norm(Ah-Ah,1) == norm(Ah-Ah,2) == norm(Ah-Ah,Inf) == 0
@test iszero(opnorm(Ah-Ah,1)) && iszero(opnorm(Ah-Ah,2)) && iszero(opnorm(Ah-Ah,Inf)) && iszero(opnorm(Ah-Ah))
@test trace(Ah-Ah) == 0 && iszero(tr(Ah-Ah))

@test Ah1i*Ah1 ≈ I ≈ Ah1*Ah1i 
@test hrchop(Ah1i; tol = 1.e-10) ≈ hrchop(Ah1i; tol = eps()) ≈ hrchop(Ah1i; tol = 1.e-10) 
@test hrtrunc(Ah1i,19) ≈ hrtrunc(Ah1i,20)

@test blockdiag(Ah,Ch)(t) ≈ bldiag(Ah(t),Ch(t))

Ah = convert(HarmonicArray,PeriodicFunctionMatrix(A,4*pi));
Ch = convert(HarmonicArray,PeriodicFunctionMatrix(C,2*pi));
Cdh = convert(HarmonicArray,PeriodicFunctionMatrix(Cd,2*pi));
Xh = convert(HarmonicArray,PeriodicFunctionMatrix(X,8*pi));
Xderh = convert(HarmonicArray,PeriodicFunctionMatrix(Xder,8*pi));
@test Ah*Xh+Xh*Ah'+Ch ≈  pmderiv(Xh) ≈ Xderh
@test Ah'*Xh+Xh*Ah+Cdh ≈ -pmderiv(Xh) 
@test norm(Ah*Xh+Xh*Ah'+Ch-pmderiv(Xh),Inf) < 1.e-7 && norm(pmderiv(Xh)- Xderh) < 1.e-7
@test norm(Ah'*Xh+Xh*Ah+Cdh+pmderiv(Xh),1) < 1.e-7

@test tpmeval(Ah,1)[1:2,1:1] == tpmeval(Ah[1:2,1],1) && lastindex(Ah,1) == 2 && lastindex(Ah,2) == 2

t = rand(); 
@test [Ah Ch](t) ≈ [Ah(t) Ch(t)]
@test [Ah; Ch](t) ≈ [Ah(t); Ch(t)]
@test blockdiag(Ah,Ch)(t) ≈ bldiag(Ah(t),Ch(t))


# PeriodicSymbolicMatrix
@show "symbolic operations"
@variables t
A11 = [0  1; -10*cos(t)-1 -24-19*sin(t)]
X1 =  [1+cos(t) 0; 0 1+sin(t)] 
X1der = [-sin(t) 0; 0 cos(t)] 
As = PeriodicSymbolicMatrix(A11,2*pi)
Cs = PeriodicSymbolicMatrix(X1der - A11*X1-X1*A11', 2*pi)
Cds = PeriodicSymbolicMatrix(-(A11'*X1+X1*A11+X1der),2*pi)
Xs = PeriodicSymbolicMatrix(X1,2*pi)
Xders = PeriodicSymbolicMatrix(X1der,2*pi)

@test set_period(set_period(As,4pi),2pi) == As
Ast = As'
@test transpose(As) == Ast
@test tr(As) == tr(Ast)
@test trace(As) == trace(Ast)
@test opnorm(As,1) == opnorm(Ast,Inf) && opnorm(As,2) == opnorm(Ast,2) && opnorm(As,Inf) == opnorm(Ast,1)
@test norm(As,1) == norm(Ast,1) && norm(As,2) == norm(Ast,2) && norm(As,Inf) == norm(Ast,Inf)
@test_throws ArgumentError opnorm(As,3)
@test_throws ArgumentError norm(As,3)

At=convert(PeriodicFunctionMatrix,As)
@test ≈(convert(PeriodicSymbolicMatrix,At),As)

Ats = convert(PeriodicTimeSeriesMatrix,At);
@test ≈(convert(PeriodicSymbolicMatrix,Ats),As)

Ah = convert(HarmonicArray,As)
@test ≈(convert(PeriodicSymbolicMatrix,Ah),As)

Ats2 = convert(PeriodicTimeSeriesMatrix,As)
@test ≈(convert(PeriodicSymbolicMatrix,Ats2),As)

At1=PeriodicFunctionMatrix(A11,2pi)
@test ≈(convert(PeriodicSymbolicMatrix,At1),As)

@test issymmetric(Cs) && issymmetric(Cds) && issymmetric(Xs) && issymmetric(Xders)
@test As*Xs+Xs*As'+Cs ==  pmderiv(Xs) == Xders
@test As'*Xs+Xs*As+Cds == -pmderiv(Xs) 
@test norm(As*Xs+Xs*As'+Cs - pmderiv(Xs),Inf) < 1.e-7 && norm(pmderiv(Xs)- Xders) < 1.e-7
@test norm(As'*Xs+Xs*As+Cds + pmderiv(Xs),1) < 1.e-7
@test As*Xs+Xs*As'+Cs ≈  pmderiv(Xs) ≈ Xders
@test As'*Xs+Xs*As+Cds ≈ -pmderiv(Xs) 
@test transpose(As) == As'
ac = rand(2,2)
@test As+ac ≈ ac+As
@test As-ac ≈ -(ac-As)
@test As-I == -(I-As)
@test As+At ≈ At+As
@test As-At ≈ -(At-As)
@test (As*ac)(1) ≈ As(1)*ac
@test (ac*As)(1) ≈ ac*As(1)
@test (As*At)(1) ≈ As(1)*At(1)
@test (At*As)(1) ≈ At(1)*As(1)
@test As/2 == 0.5*As
@test As*I == I*As
@test [As ac](1) == [As(1) ac] && [ac As](1) == [ac As(1)]
@test [As; ac](1) == [As(1); ac] && [ac; As](1) == [ac; As(1)]
@test horzcat(As,ac)(1) == [As(1) ac] && horzcat(ac,As)(1) == [ac As(1)]
@test vertcat(As,ac)(1) == [As(1); ac] && vertcat(ac,As)(1) == [ac; As(1)]


D = rand(2,2)
@test As+I == I+As && As*5 == 5*As && As*D ≈ -As*(-D) && iszero(As-As) && !iszero(As)
@test PeriodicSymbolicMatrix(D,2*pi) == PeriodicSymbolicMatrix(D,4*pi,nperiod=2) && 
      PeriodicSymbolicMatrix(D,2*pi) ≈ PeriodicSymbolicMatrix(D,4*pi,nperiod=2)

As = PeriodicSymbolicMatrix(A11,4*pi,nperiod=2)
Cs = PeriodicSymbolicMatrix(X1der - A11*X1-X1*A11', 2*pi)
Cds = PeriodicSymbolicMatrix(-(A11'*X1+X1*A11+X1der),2*pi)
Xs = PeriodicSymbolicMatrix(X1,8*pi,nperiod=4)
Xders = PeriodicSymbolicMatrix(X1der,16*pi,nperiod = 8)
@test As*Xs+Xs*As'+Cs ==  pmderiv(Xs) == Xders
@test As'*Xs+Xs*As+Cds == -pmderiv(Xs)
@test As*Xs+Xs*As'+Cs ≈  pmderiv(Xs) ≈ Xders
@test As'*Xs+Xs*As+Cds ≈ -pmderiv(Xs)

@test (Symbolics.simplify(inv(As)*As) ≈ I) && (I ≈ Symbolics.simplify(As*inv(As))) 
@test inv(As)*As == I == As*inv(As) 

@test iszero(As[1:2,1:1].F - As.F[1:2,1:1]) && lastindex(As,1) == size(As,1) && lastindex(As,2) == size(As,2)

t = rand(); 
@test [As Cs](t) ≈ [As(t) Cs(t)]
@test [As; Cs](t) ≈ [As(t); Cs(t)]
@test blockdiag(As,Cs)(t) ≈ bldiag(As(t),Cs(t))

# FourierFunctionMatrix
@show "Fourier series operations"
@time Af = convert(FourierFunctionMatrix,PeriodicFunctionMatrix(A,2*pi));
@time Af1 = convert(FourierFunctionMatrix,PeriodicFunctionMatrix(A1,2*pi));
@time Af2 = convert(FourierFunctionMatrix,PeriodicFunctionMatrix(A,4*pi));
@time Cf = convert(FourierFunctionMatrix,PeriodicFunctionMatrix(C,2*pi));
@time Cdf = convert(FourierFunctionMatrix,PeriodicFunctionMatrix(Cd,2*pi));
@time Xf = convert(FourierFunctionMatrix,PeriodicFunctionMatrix(X,2*pi));
@time Xderf = convert(FourierFunctionMatrix,PeriodicFunctionMatrix(Xder,2*pi));

@test set_period(set_period(Af,4pi),2pi) == Af
Aft = Af'
@test transpose(Af) == Aft
@test tr(Af) == tr(Aft)
@test trace(Af) == trace(Aft)
@test opnorm(Af,2) ≈ opnorm(Aft,2) 
@test norm(Af,1) == norm(Aft,1) && norm(Af,2) == norm(Aft,2) && norm(Af,Inf) == norm(Aft,Inf)
@test_throws ArgumentError opnorm(Af,1)
@test_throws ArgumentError norm(Af,3)

@test FourierFunctionMatrix(Af.M) == Af
a = rand(2,2);
Af11=FourierFunctionMatrix(a,3)
@test a == pmaverage(Af11)
@test isconstant(FourierFunctionMatrix(rand(1,1),3))

As = PeriodicSymbolicMatrix(A11,2*pi)

@time At1=convert(PeriodicFunctionMatrix,Af) 
@test ≈(convert(FourierFunctionMatrix,At1),Af)

@time Ah1 = convert(HarmonicArray,Af)
@test ≈(convert(FourierFunctionMatrix,Ah1),Af)

@time Ah1 = convert(HarmonicArray,Af*Af)
@test ≈(convert(FourierFunctionMatrix,Ah1),Af*Af)

@time Ats1 = convert(PeriodicTimeSeriesMatrix,Af);
@test ≈(convert(FourierFunctionMatrix,Ats1),Af)

@time Ats2 = convert(PeriodicTimeSeriesMatrix,Af)
@test ≈(convert(FourierFunctionMatrix,Ats2),Af)

@time Af11 = convert(FourierFunctionMatrix,As)
@time As1 = convert(PeriodicSymbolicMatrix,Af11)

@test ≈(As1,As) 
@test ≈(Af11,Af)
@test ≈(As,PeriodicSymbolicMatrix(ffm2psm(Af,0:0)+ffm2psm(Af,1:1),As.period))
@test domain((Af+Af2).M).b ≈ 4pi && domain((Af*Af2).M).b ≈ 4pi

ac = rand(2,2)
@test Af+ac ≈ ac+Af
@test Af-ac ≈ -(ac-Af)
@test Af-I == -(I-Af)
@test Af+At ≈ At+Af
@test Af-At ≈ -(At-Af)
@test (Af*ac)(1) ≈ Af(1)*ac
@test (ac*Af)(1) ≈ ac*Af(1)
@test (Af*At)(1) ≈ Af(1)*At(1)
@test (At*Af)(1) ≈ At(1)*Af(1)
@test Af/2 == 0.5*Af
@test Af*I == I*Af
@test [Af ac](1) ≈ [Af(1) ac] && [ac Af](1) ≈ [ac Af(1)]
@test [Af; ac](1) ≈ [Af(1); ac] && [ac; Af](1) ≈ [ac; Af(1)]
@test horzcat(Af,ac)(1) ≈ [Af(1) ac] && horzcat(ac,Af)(1) ≈ [ac Af(1)]
@test vertcat(Af,ac)(1) ≈ [Af(1); ac] && vertcat(ac,Af)(1) ≈ [ac; Af(1)]


@test issymmetric(Cf) && issymmetric(Cdf) && issymmetric(Xf) && issymmetric(Xderf)
@test Af*Xf+Xf*Af'+Cf ≈  pmderiv(Xf) ≈ Xderf
@test Af'*Xf+Xf*Af+Cdf ≈ -pmderiv(Xf) 
@test norm(Af*Xf+Xf*Af'+Cf-pmderiv(Xf),Inf) < 1.e-7 && norm(pmderiv(Xf)- Xderf) < 1.e-7
@test norm(Af'*Xf+Xf*Af+Cdf+pmderiv(Xf),1) < 1.e-7

D = rand(2,2)
@test Af+I == I+Af && Af*5 == 5*Af && Af*D ≈ -Af*(-D)  && iszero(Af-Af) && !iszero(Af) 
@test FourierFunctionMatrix(D,2*pi) == FourierFunctionMatrix(D,4*pi) && 
      FourierFunctionMatrix(D,2*pi) ≈ FourierFunctionMatrix(D,4*pi)
@test inv(Af1)*Af1 ≈ I ≈ Af1*inv(Af1) 
@test norm(Af-Af,1) == norm(Af-Af,2) == norm(Af-Af,Inf) == 0
@test iszero(opnorm(Af-Af,2)) && iszero(opnorm(Af-Af)) 
@test trace(Af-Af) == 0 && iszero(tr(Af-Af))

@test blockdiag(Af,Cf)(t) ≈ bldiag(Af(t),Cf(t))

ts = sort(rand(10))
@test  tvmeval(Af,ts) ≈ tpmeval.(Ref(Af),ts)


@time Af = convert(FourierFunctionMatrix,PeriodicFunctionMatrix(A,4*pi,nperiod=2));
@time Cf = convert(FourierFunctionMatrix,PeriodicFunctionMatrix(C,2*pi));
@time Cdf = convert(FourierFunctionMatrix,PeriodicFunctionMatrix(Cd,2*pi));
@time Xf = convert(FourierFunctionMatrix,PeriodicFunctionMatrix(X,8*pi,nperiod=4));
@time Xderf = convert(FourierFunctionMatrix,PeriodicFunctionMatrix(Xder,16*pi,nperiod=8));
@test Af*Xf+Xf*Af'+Cf ≈  pmderiv(Xf) ≈ Xderf
@test Af'*Xf+Xf*Af+Cdf ≈ -pmderiv(Xf) 
@test norm(Af*Xf+Xf*Af'+Cf-pmderiv(Xf),Inf) < 1.e-7 && norm(pmderiv(Xf)- Xderf) < 1.e-7
@test norm(Af'*Xf+Xf*Af+Cdf+pmderiv(Xf),1) < 1.e-7

@test tpmeval(Af,1)[1:2,1:1] == tpmeval(Af[1:2,1],1) && lastindex(Af,1) == 2 && lastindex(Af,2) == 2

t = rand(); 
@test [Af Cf](t) ≈ [Af(t) Cf(t)]
@test [Af; Cf](t) ≈ [Af(t); Cf(t)]
@test blockdiag(Af,Cf)(t) ≈ bldiag(Af(t),Cf(t))

# PeriodicTimeSeriesMatrix
t1 = -rand(Float32,2,2); t2 = rand(Float32,2,2); Ats1 = PeriodicTimeSeriesMatrix{:c,Float64}([t1,t2],2)
@test Ats1.ts == [0.,1.]
Ats = convert(PeriodicTimeSeriesMatrix,PeriodicFunctionMatrix(A,2*pi));
Cts = convert(PeriodicTimeSeriesMatrix,PeriodicFunctionMatrix(C,2*pi));
Cdts = convert(PeriodicTimeSeriesMatrix,PeriodicFunctionMatrix(Cd,2*pi));
Xts = convert(PeriodicTimeSeriesMatrix,PeriodicFunctionMatrix(X,2*pi));
Xderts = convert(PeriodicTimeSeriesMatrix,PeriodicFunctionMatrix(Xder,2*pi));
@test issymmetric(Cts) && issymmetric(Cdts) && issymmetric(Xts) && issymmetric(Xderts)

@test set_period(set_period(Ats,4pi),2pi) == Ats
@test Ats[:,:] == Ats
@test lastindex(Ats,1) == size(Ats,1) && lastindex(Ats,2) == size(Ats,2)
@test Ats[1] == Ats(0) == Ats[end+1]


# @test Ats*Xts+Xts*Ats'+Cts ≈  pmderiv(Xts) ≈ Xderts
# @test Ats'*Xts+Xts*Ats+Cdts ≈ -pmderiv(Xts) 
@test norm(Ats*Xts+Xts*Ats'+Cts-pmderiv(Xts),Inf) < 1.e-6 && norm(Ats*Xts+Xts*Ats'+Cts-Xderts,Inf) < 1.e-7
@test norm(Ats'*Xts+Xts*Ats+Cdts+pmderiv(Xts),Inf) < 1.e-6 && norm(Ats'*Xts+Xts*Ats+Cdts+Xderts,Inf) < 1.e-6

D = rand(2,2)
@test Ats+I == I+Ats && Ats*5 == 5*Ats && Ats*D ≈ -Ats*(-D)  && iszero(Ats-Ats) && !iszero(Ats)
@test PeriodicTimeSeriesMatrix(D,2*pi) !== PeriodicTimeSeriesMatrix(D,4*pi) && 
      !(PeriodicTimeSeriesMatrix(D,2*pi) ≈ PeriodicTimeSeriesMatrix(D,4*pi))
@test inv(Ats)*Ats ≈ I ≈ Ats*inv(Ats) 
@test norm(Ats-Ats,1) == norm(Ats-Ats,2) == norm(Ats-Ats,Inf) == 0
@test iszero(opnorm(Ats-Ats,1)) && iszero(opnorm(Ats-Ats,2)) && iszero(opnorm(Ats-Ats,Inf)) && 
      iszero(opnorm(Ats-Ats)) 
@test trace(Ats-Ats) == 0 && iszero(tr(Ats-Ats))

t = rand(); 
@test blockdiag(Ats,Cts)(t) ≈ bldiag(Ats(t),Cts(t))

# same tsub
# TA = PeriodicTimeSeriesMatrix([[i;;] for i in 1:4],30,nperiod=6)
# TB = PeriodicTimeSeriesMatrix([[i;;] for i in 1:2],10,nperiod=2)
TA = PeriodicTimeSeriesMatrix([rand(2,2) for i in 1:4],30,nperiod=6)
TB = PeriodicTimeSeriesMatrix([rand(2,2) for i in 1:2],10,nperiod=2)
AB = TA+TB
ts = sort(rand(100)*AB.period); 
@test norm(AB.(ts).-TA.(ts).-TB.(ts)) < 1.e-10
AB = TA*TB
ts = sort(rand(10)*AB.period); 
@test norm(AB.(ts).-TA.(ts).*TB.(ts)) < 1.e-10


# same period
# TA = PeriodicTimeSeriesMatrix([[i;;] for i in 1:4],30,nperiod=6)
# TB = PeriodicTimeSeriesMatrix([[i;;] for i in 1:2],30,nperiod=2)
TA = PeriodicTimeSeriesMatrix([rand(2,2) for i in 1:4],30,nperiod=6)
TB = PeriodicTimeSeriesMatrix([rand(2,2) for i in 1:2],30,nperiod=2)
AB = TA+TB
ts = sort(rand(10)*AB.period); 
@test norm(AB.(ts).-TA.(ts).-TB.(ts)) < 1.e-10
AB = TA*TB
@test norm(AB.(ts).-TA.(ts).*TB.(ts)) < 1.e-10



Ats = convert(PeriodicTimeSeriesMatrix,PeriodicFunctionMatrix(A,4*pi,nperiod=2));
Cts = convert(PeriodicTimeSeriesMatrix,PeriodicFunctionMatrix(C,2*pi));
Cdts = convert(PeriodicTimeSeriesMatrix,PeriodicFunctionMatrix(Cd,2*pi));
Xts = convert(PeriodicTimeSeriesMatrix,PeriodicFunctionMatrix(X,8*pi,nperiod=4));
Xderts = convert(PeriodicTimeSeriesMatrix,PeriodicFunctionMatrix(Xder,16*pi,nperiod=8));
@test norm(Ats*Xts+Xts*Ats'+Cts-pmderiv(Xts),Inf) < 1.e-6 && norm(Ats*Xts+Xts*Ats'+Cts-Xderts,Inf) < 1.e-6 
@test norm(Ats'*Xts+Xts*Ats+Cdts+pmderiv(Xts),Inf) < 1.e-6 && norm(Ats'*Xts+Xts*Ats+Cdts+Xderts,Inf) < 1.e-6

@test Ats[1:2,1:1].values == [Ats.values[i][1:2,1:1] for i in 1:length(Ats)] && lastindex(Ats,1) == 2 && lastindex(Ats,2) == 2

t = rand(); 
@test [Ats Cts](t) ≈ [Ats(t) Cts(t)]
@test [Ats; Cts](t) ≈ [Ats(t); Cts(t)]
@test blockdiag(Ats,Cts)(t) ≈ bldiag(Ats(t),Cts(t))


# PeriodicSwitchingMatrix
t1 = -rand(2,2); t2 = rand(2,2); Asw = PeriodicSwitchingMatrix([t1,t2],[0.,1.],2)
t1 = -rand(Float32,2,2); t2 = rand(Float32,2,2); Asw1 = PeriodicSwitchingMatrix{:c,Float64}([t1,t2],[0.,1.],2)
Asw2 = PeriodicSwitchingMatrix(rand(2,2,2),[0.,1.],2)

@test set_period(set_period(Asw,4),2) == Asw
@test Asw[:,:] == Asw
@test lastindex(Asw,1) == size(Asw,1) && lastindex(Asw,2) == size(Asw,2)

@test convert(PeriodicFunctionMatrix,Asw)(1) ≈ convert(PeriodicFunctionMatrix{:c,BigFloat},Asw)(1)
@test convert(HarmonicArray,Asw)(1) ≈ Asw(1)


t1 = rand(2,2); t2 = rand(2,2); Csw = PeriodicSwitchingMatrix([t1,t2],[0.,1.5],2)
@test Csw == convert(PeriodicSwitchingMatrix,convert(PeriodicTimeSeriesMatrix,Csw,ns=10))
t1 = rand(2,2); t2 = rand(2,2); Tsw = PeriodicSwitchingMatrix([t1,t2],[0.,pi/2],2)
@test norm((Tsw - convert(PeriodicSwitchingMatrix,convert(PeriodicTimeSeriesMatrix,Tsw,ns=10))).(rand(10*2))) == 0

t = 2*rand(); 
@test (Asw+Csw)(t) ≈ Asw(t)+Csw(t)
@test (Asw*Csw)(t) ≈ Asw(t)*Csw(t)
@test [Asw Csw](t) ≈ [Asw(t) Csw(t)]
@test [Asw; Csw](t) ≈ [Asw(t); Csw(t)]
@test norm(Asw-2*Asw+Asw) == 0
D = rand(2,2)
@test Asw+I == I+Asw && Asw*5 == 5*Asw && Asw*D ≈ -Asw*(-D)  && iszero(Asw-Asw) && !iszero(Asw)
@test PeriodicSwitchingMatrix(D,2) !== PeriodicSwitchingMatrix(D,4) && 
      !(PeriodicSwitchingMatrix(D,2) ≈ PeriodicSwitchingMatrix(D,4))
@test issymmetric(Csw*Csw')
@test inv(Asw)*Asw ≈ I ≈ Asw*inv(Asw) 
@test norm(Asw-Asw,1) == norm(Asw-Asw,2) == norm(Asw-Asw,Inf) == 0
@test iszero(opnorm(Asw-Asw,1)) && iszero(opnorm(Asw-Asw,2)) && iszero(opnorm(Asw-Asw,Inf)) && 
      iszero(opnorm(Asw-Asw)) 
@test trace(Asw-Asw) == 0 && iszero(tr(Asw-Asw))
t = rand(); 
@test blockdiag(Asw,Csw)(t) ≈ bldiag(Asw(t),Csw(t))




# PeriodicArray
n = 5; pa = 3; px = 6;   
Ad = 0.5*PeriodicArray(rand(Float64,n,n,pa),pa);
Ad1 = PeriodicArray(Ad.M,pa);
@test Ad(0) == Ad.M[:,:,1]
@test getpm(Ad,4) == getpm(Ad,1)
x = rand(n,n,px); [x[:,:,i] = x[:,:,i]'+x[:,:,i] for i in 1:px];
Xd = PeriodicArray(x,px);
Qdf = -Ad*Xd*Ad'+pmshift(Xd); Qdf = (Qdf+transpose(Qdf))/2
Qdr = -Ad'*pmshift(Xd)*Ad+Xd; Qdr = (Qdr+transpose(Qdr))/2

@test set_period(set_period(Ad,2*pa),pa) == Ad
@test propertynames(Ad) == (:dperiod, :Ts, :M, :period, :nperiod)
@test lastindex(Ad) == Ad.dperiod

@test !iscontinuous(Ad) 

@test convert(PeriodicTimeSeriesMatrix,Ad) == convert(PeriodicTimeSeriesMatrix,convert(PeriodicMatrix,Ad))
@test convert(PeriodicArray,convert(SwitchingPeriodicMatrix,Ad)) == Ad

@test reverse(reverse(Ad)) == Ad



# Xf = pfdlyap(Ad, Qdf);
@test Ad*Xd*Ad' + Qdf ≈ pmshift(Xd) 
# Xr = prdlyap(Ad, Qdr);
@test Ad'*pmshift(Xd)*Ad + Qdr ≈ Xd 


Qds = pmshift(Qdf); 

@test issymmetric(Qdf) && issymmetric(Qds) && isequal(pmshift(pmshift(Qdf,1),-1),Qdf) && iszero(Qdf-Qdf')
@test inv(Ad)*Ad ≈ I ≈ Ad*inv(Ad) && Ad+I == I+Ad
@test norm(Ad-Ad,1) == norm(Ad-Ad,2) == norm(Ad-Ad,Inf) == 0
@test_throws ArgumentError norm(Ad,3)
@test iszero(opnorm(Ad-Ad,1)) && iszero(opnorm(Ad-Ad,2)) && iszero(opnorm(Ad-Ad,Inf)) && iszero(opnorm(Ad-Ad))
@test trace(Ad-Ad) == 0 && iszero(tr(Ad-Ad))
aa = rand(5,5)
@test Ad+aa ≈ aa+Ad
@test Ad-aa ≈ -(aa-Ad)
@test Ad-I ≈ -(I-Ad)
@test (Ad+Ad')/2 ≈ pmsymadd!(Ad1,0.5) && (Ad+Ad') ≈ pmsymadd!(Ad1) && issymmetric(Ad1)
@test (Ad*aa)(1) ≈ Ad(1)*aa && (aa*Ad)(1) ≈ aa*Ad(1)
@test Ad*I == I*Ad
Adt = transpose(Ad)
@test pmmuladdsym(Ad1, Ad, Adt, 1, 1) ≈ Ad1+Ad*Adt
@test pmmultraddsym(Ad1, Ad, Ad, 1, 1) ≈ Ad1+Adt*Ad
@test pmmuladdtrsym(Ad1, Ad, Ad, 1, 1) ≈ Ad1+Ad*Adt
@test pmmulsym(Ad, Adt, 1) ≈ Ad*Adt
@test pmmultrsym(Ad, Ad, 1) ≈ Adt*Ad
@test pmmultrsym(Ad, Ad, 1) ≈ Ad*Adt
@test [Ad aa](1) ≈ [Ad(1) aa] && [aa Ad](1) ≈ [aa Ad(1)]
@test [Ad; aa](1) ≈ [Ad(1); aa] && [aa; Ad](1) ≈ [aa; Ad(1)]
@test horzcat(Ad,aa)(1) ≈ [Ad(1) aa] && horzcat(aa,Ad)(1) ≈ [aa Ad(1)]
@test vertcat(Ad,aa)(1) ≈ [Ad(1); aa] && vertcat(aa,Ad)(1) ≈ [aa; Ad(1)]


D = rand(n,n)
@test Ad*5 == 5*Ad && Ad*D ≈ -Ad*(-D)  && iszero(Ad-Ad) && !iszero(Ad)
@test PeriodicArray(D,2*pi) !== PeriodicArray(D,4*pi) && 
      !(PeriodicArray(D,2*pi) ≈ PeriodicArray(D,4*pi))

@test blockdiag(Ad,Xd)[10] ≈ bldiag(Ad[10],Xd[10])      


Ad1 = PeriodicArray(Ad.M,2*pa;nperiod=2);
Ad2 = set_period(Ad,2*pa)
Xd1 = PeriodicArray(x,3*px; nperiod = 3);
Qdf1 = -Ad1*Xd1*Ad1'+pmshift(Xd1); Qdf1 = (Qdf1+transpose(Qdf1))/2
Qdr1 = -Ad1'*pmshift(Xd1)*Ad1+Xd1; Qdr1 = (Qdr1+transpose(Qdr1))/2

@test Ad1 == Ad2 && (Ad1+Ad)/3 ≈ 2*Ad2/3 && issymmetric(Qdf1+Qdr1)
# Xf1 = pfdlyap(Ad1, Qdf1);
@test Ad*Xd*Ad' + Qdf ≈ pmshift(Xd) 
# Xr1 = prdlyap(Ad1, Qdr1);
@test Ad'*pmshift(Xd)*Ad + Qdr ≈ Xd 

@test Ad[1:2,1:1].M == Ad.M[1:2,1:1,:]  && lastindex(Ad,1) == n && lastindex(Ad,2) == n

@test blockdiag(Ad1,Xd1)[10] ≈ bldiag(Ad1[10],Xd1[10])  


# PeriodicMatrix
n = 5; pa = 3; px = 6;   
Ad = 0.5*PeriodicMatrix([rand(Float64,n,n) for i in 1:pa],pa);
Ad1 = PeriodicMatrix(Ad.M,pa);
@test Ad(0) == Ad.M[1]
@test getpm(Ad,4) == getpm(Ad,1)
x = [rand(n,n) for i in 1:px]
Xd = PeriodicMatrix([ x[i]+x[i]' for i in 1:px],px);
Qdf = -Ad*Xd*Ad'+pmshift(Xd); Qdf = (Qdf+transpose(Qdf))/2
Qdr = -Ad'*pmshift(Xd)*Ad+Xd; Qdr = (Qdr+transpose(Qdr))/2

@test set_period(set_period(Ad,2*pa),pa) == Ad
@test propertynames(Ad) == (:dperiod, :Ts, :M, :period, :nperiod)
@test lastindex(Ad) == Ad.dperiod

#Xf = pfdlyap(Ad, Qdf);
@test Ad*Xd*Ad' + Qdf ≈ pmshift(Xd) 
# Xr = prdlyap(Ad, Qdr);
@test Ad'*pmshift(Xd)*Ad + Qdr ≈ Xd 


Qds = pmshift(Qdf); 

@test issymmetric(Qdf) && issymmetric(Qds) && isequal(pmshift(pmshift(Qdf,1),-1),Qdf) && iszero(Qdf-Qdf')
@test inv(Ad)*Ad ≈ I ≈ Ad*inv(Ad) && Ad+I == I+Ad
@test norm(Ad-Ad,1) == norm(Ad-Ad,2) == norm(Ad-Ad,Inf) == 0
@test iszero(opnorm(Ad-Ad,1)) && iszero(opnorm(Ad-Ad,2)) && iszero(opnorm(Ad-Ad,Inf))
@test trace(Ad-Ad) == 0 && iszero(tr(Ad-Ad))
@test reverse(reverse(Ad)) == Ad
@test norm(Ad-Ad,1) == norm(Ad-Ad,2) == norm(Ad-Ad,Inf) == 0
@test_throws ArgumentError norm(Ad,3)
@test iszero(opnorm(Ad-Ad,1)) && iszero(opnorm(Ad-Ad,2)) && iszero(opnorm(Ad-Ad,Inf)) && iszero(opnorm(Ad-Ad))
aa = rand(5,5)
@test Ad+aa ≈ aa+Ad
@test Ad-aa ≈ -(aa-Ad)
@test Ad-I ≈ -(I-Ad)
@test (Ad+Ad')/2 ≈ pmsymadd!(Ad1,0.5) && (Ad+Ad') ≈ pmsymadd!(Ad1) && issymmetric(Ad1)
@test (Ad*aa)(1) ≈ Ad(1)*aa && (aa*Ad)(1) ≈ aa*Ad(1)
@test Ad*I == I*Ad
Adt = transpose(Ad)
@test pmmuladdsym(Ad1, Ad, Adt, 1, 1) ≈ Ad1+Ad*Adt
@test pmmultraddsym(Ad1, Ad, Ad, 1, 1) ≈ Ad1+Adt*Ad
@test pmmuladdtrsym(Ad1, Ad, Ad, 1, 1) ≈ Ad1+Ad*Adt
@test pmmulsym(Ad, Adt, 1) ≈ Ad*Adt
@test pmmultrsym(Ad, Ad, 1) ≈ Adt*Ad
@test pmmultrsym(Ad, Ad, 1) ≈ Ad*Adt
@test [Ad aa](1) ≈ [Ad(1) aa] && [aa Ad](1) ≈ [aa Ad(1)]
@test [Ad; aa](1) ≈ [Ad(1); aa] && [aa; Ad](1) ≈ [aa; Ad(1)]
@test horzcat(Ad,aa)(1) ≈ [Ad(1) aa] && horzcat(aa,Ad)(1) ≈ [aa Ad(1)]
@test vertcat(Ad,aa)(1) ≈ [Ad(1); aa] && vertcat(aa,Ad)(1) ≈ [aa; Ad(1)]



@test convert(PeriodicArray,Ad) == convert(PeriodicArray{:d,Float64},Ad)
@test convert(PeriodicTimeSeriesMatrix,Ad) == convert(PeriodicTimeSeriesMatrix,convert(PeriodicArray,Ad))


D = rand(n,n)
@test Ad*5 == 5*Ad  && Ad*D ≈ -Ad*(-D) && iszero(Ad-Ad) && !iszero(Ad)
@test PeriodicMatrix(D,2*pi) !== PeriodicMatrix(D,4*pi) && 
      !(PeriodicMatrix(D,2*pi) ≈ PeriodicMatrix(D,4*pi))
     
@test blockdiag(Ad,Xd)[10] ≈ bldiag(Ad[10],Xd[10])      
      


Ad1 = PeriodicMatrix(Ad.M,2*pa;nperiod=2);
Xd1 = PeriodicMatrix(Xd.M,3*px; nperiod = 3);
Qdf1 = -Ad1*Xd1*Ad1'+pmshift(Xd1); Qdf1 = (Qdf1+transpose(Qdf1))/2
Qdr1 = -Ad1'*pmshift(Xd1)*Ad1+Xd1; Qdr1 = (Qdr1+transpose(Qdr1))/2

@test Ad1 !== Ad && (Ad1+Ad)/3 ≈ 2*Ad1/3 && issymmetric(Qdf1+Qdr1)
# Xf1 = pfdlyap(Ad1, Qdf1);
@test Ad*Xd1*Ad' + Qdf ≈ pmshift(Xd1) 
#Xr1 = prdlyap(Ad1, Qdr1);
@test Ad'*pmshift(Xd1)*Ad + Qdr ≈ Xd1 

@test Ad[1:2,1:1].M == [Ad.M[i][1:2,1:1] for i in 1:length(Ad)] && lastindex(Ad,1) == [n;n;n] && lastindex(Ad,2) == [n;n;n]


@test [[Ad Ad]; [Ad Ad]] == [[Ad;Ad] [Ad;Ad]]
@test blockdiag(Ad1,Xd1)[10] ≈ bldiag(Ad1[10],Xd1[10]) 


# time-varying dimensions
na = [5, 3, 3, 4, 1]; ma = [3, 3, 4, 1, 5]; pa = 5; px = 5;   
#na = 5*na; ma = 5*ma;
Ad = PeriodicMatrix([rand(Float64,ma[i],na[i]) for i in 1:pa],pa);  
x = [rand(na[i],na[i]) for i in 1:px]
Xd = PeriodicMatrix([ x[i]+x[i]' for i in 1:px],px);
Qdf = -Ad*Xd*Ad'+pmshift(Xd); Qdf = (Qdf+transpose(Qdf))/2
Qdr = -Ad'*pmshift(Xd)*Ad+Xd; Qdr = (Qdr+transpose(Qdr))/2

@test Ad.M + pmzeros(ma,na) == Ad.M


# Xf = pfdlyap(Ad, Qdf);
@test Ad*Xd*Ad' + Qdf ≈ pmshift(Xd) 
# Xr = prdlyap(Ad, Qdr);
@test Ad'*pmshift(Xd)*Ad + Qdr ≈ Xd 

Qds = pmshift(Qdf); 

@test issymmetric(Qdf) && issymmetric(Qds) && isequal(pmshift(pmshift(Qdf,1),-1),Qdf) && iszero(Qdf-Qdf')

@test Ad*5 == 5*Ad &&  iszero(Ad-Ad) && !iszero(Ad) && Qdf + I == I+Qdf
@test_throws DimensionMismatch Ad ≈ I && Ad-I 

@test blockdiag(Ad,Xd)[10] ≈ bldiag(Ad[10],Xd[10])      



Ad1 = PeriodicMatrix(Ad.M,2*pa;nperiod=2);
Xd1 = PeriodicMatrix(Xd.M,3*px; nperiod = 3);
Qdf1 = -Ad1*Xd1*Ad1'+pmshift(Xd1); Qdf1 = (Qdf1+transpose(Qdf1))/2
Qdr1 = -Ad1'*pmshift(Xd1)*Ad1+Xd1; Qdr1 = (Qdr1+transpose(Qdr1))/2

@test Ad1 !== Ad && (Ad1+Ad)/3 ≈ 2*Ad1/3 && issymmetric(Qdf1)&& issymmetric(Qdr1)
# Xf1 = pfdlyap(Ad1, Qdf1);
@test Ad*Xd1*Ad' + Qdf ≈ pmshift(Xd1) 
# Xr1 = prdlyap(Ad1, Qdr1);
@test Ad'*pmshift(Xd1)*Ad + Qdr ≈ Xd1 

@test Ad[1:1,1:1].M == [Ad.M[i][1:1,1:1] for i in 1:length(Ad)] && lastindex(Ad,1) == size(Ad,1) && lastindex(Ad,2) == size(Ad,2)

@test blockdiag(Ad1,Xd1)[10] ≈ bldiag(Ad1[10],Xd1[10]) 


# SwitchingPeriodicMatrix
n = 2; pa = 3; px = 6; T = 10; 
Ad = 0.5*SwitchingPeriodicMatrix([rand(Float64,n,n) for i in 1:pa],[10,15,20],T);
@test Ad(0) == Ad.M[1]

x = [rand(n,n) for i in 1:px]
Xd = SwitchingPeriodicMatrix([ x[i]+x[i]' for i in 1:px],[2, 3, 5,7, 9, 10],T;nperiod=2);
@test Ad.Ts == Xd.Ts

@test set_period(set_period(Ad,2*T),T) == Ad
@test propertynames(Ad) == (:dperiod, :Ts, :M, :ns, :period, :nperiod)
@test lastindex(Ad) == Ad.dperiod


@test Ad == pmshift(pmshift(Ad),-1)
@test Ad == pmshift(pmshift(Ad,10),-10)

Qdf = -Ad*Xd*Ad'+pmshift(Xd); Qdf = (Qdf+transpose(Qdf))/2
Qdr = -Ad'*pmshift(Xd)*Ad+Xd; Qdr = (Qdr+transpose(Qdr))/2

# Xf = pfdlyap(Ad, Qdf);
@test norm(Ad*Xd*Ad' + Qdf - pmshift(Xd)) < 1.e-7 
# Xr = prdlyap(Ad, Qdr);
@test norm(Ad'*pmshift(Xd)*Ad + Qdr - Xd) < 1.e-7 



@test issymmetric(Xd) && iszero(Xd-Xd')
@test inv(Ad)*Ad ≈ I ≈ Ad*inv(Ad) && Ad+I == I+Ad
@test Ad == reverse(reverse(Ad))
@test Ad == convert(SwitchingPeriodicMatrix,convert(PeriodicMatrix,Ad))
@test Ad ≈ convert(SwitchingPeriodicMatrix,convert(PeriodicMatrix,Ad))
@test norm(Ad,Inf) == norm(convert(PeriodicMatrix,Ad),Inf)
@test norm(Ad,1) ≈ norm(convert(PeriodicMatrix,Ad),1)
@test norm(Ad,2) ≈ norm(convert(PeriodicMatrix,Ad),2)
@test iszero(opnorm(Ad-Ad,Inf)) && iszero(opnorm(Ad-Ad,1)) && iszero(opnorm(Ad-Ad,2))
@test trace(Ad) ≈ trace(convert(PeriodicMatrix,Ad)) && tr(Ad) ≈ convert(SwitchingPeriodicMatrix,tr(convert(PeriodicMatrix,Ad)))

@test convert(PeriodicArray,Ad) == convert(PeriodicArray,convert(PeriodicMatrix,Ad))

D = rand(n,n)
@test Ad*5 == 5*Ad  && Ad*D ≈ -Ad*(-D) && iszero(Ad-Ad) && !iszero(Ad)
@test SwitchingPeriodicMatrix(D,2*pi) !== SwitchingPeriodicMatrix(D,4*pi) && 
      !(SwitchingPeriodicMatrix(D,2*pi) ≈ SwitchingPeriodicMatrix(D,4*pi))


@test Ad[1:2,1:1].M == [Ad.M[i][1:2,1:1] for i in 1:length(Ad.M)] && lastindex(Ad,1) == size(Ad,1) && lastindex(Ad,2) == size(Ad,2)

@test [[Ad Xd]; [Xd Ad]] == [[Ad;Xd] [Xd;Ad]]

@test blockdiag(Ad,Xd)[10] ≈ bldiag(Ad[10],Xd[10])   

# SwitchingPeriodicArray
n = 2; pa = 3; px = 6; T = 10; 
Ad = 0.5*SwitchingPeriodicArray(rand(Float64,n,n,pa),[10,15,20],T);
@test Ad(0) == Ad.M[:,:,1]

x = pmsymadd!(PeriodicArray(rand(n,n,px),T));
Xd = SwitchingPeriodicArray(x.M,[2, 3, 5,7, 9, 10],T;nperiod=2);
@test Ad.Ts == Xd.Ts

@test set_period(set_period(Ad,2*T),T) == Ad
@test propertynames(Ad) == (:dperiod, :Ts, :M, :ns, :period, :nperiod)
@test lastindex(Ad) == Ad.dperiod
@test size(Ad) == (size(Ad,1),size(Ad,2))


@test Ad == pmshift(pmshift(Ad),-1)
@test Ad == pmshift(pmshift(Ad,10),-10)

Qdf = -Ad*Xd*Ad'+pmshift(Xd); Qdf = (Qdf+transpose(Qdf))/2
Qdr = -Ad'*pmshift(Xd)*Ad+Xd; Qdr = (Qdr+transpose(Qdr))/2

#Xf = pfdlyap(Ad, Qdf);
@test norm(Ad*Xd*Ad' + Qdf - pmshift(Xd)) < 1.e-7
# Xr = prdlyap(Ad, Qdr);
@test norm(Ad'*pmshift(Xd)*Ad + Qdr - Xd) < 1.e-7 



@test issymmetric(Xd) && iszero(Xd-Xd')
@test inv(Ad)*Ad ≈ I ≈ Ad*inv(Ad) && Ad+I == I+Ad
@test Ad == reverse(reverse(Ad))
@test Ad == convert(SwitchingPeriodicArray,convert(PeriodicArray,Ad))
@test Ad ≈ convert(SwitchingPeriodicArray,convert(PeriodicArray,Ad))
@test norm(Ad,Inf) == norm(convert(PeriodicArray,Ad),Inf)
@test norm(Ad,1) ≈ norm(convert(PeriodicArray,Ad),1)
@test norm(Ad,2) ≈ norm(convert(PeriodicArray,Ad),2)
@test iszero(opnorm(Ad-Ad,Inf)) && iszero(opnorm(Ad-Ad,1)) && iszero(opnorm(Ad-Ad,2))
@test trace(Ad) ≈ trace(convert(PeriodicArray,Ad)) && tr(Ad) ≈ convert(SwitchingPeriodicArray,tr(convert(PeriodicArray,Ad)))

@test convert(PeriodicMatrix,convert(PeriodicArray,Ad)) == convert(PeriodicMatrix,Ad)
@test convert(SwitchingPeriodicArray,convert(PeriodicArray,Ad)) == Ad
@test convert(SwitchingPeriodicArray,convert(SwitchingPeriodicMatrix,Ad)) == Ad


D = rand(n,n)
@test Ad*5 == 5*Ad  && Ad*D ≈ -Ad*(-D) && iszero(Ad-Ad) && !iszero(Ad)
@test SwitchingPeriodicArray(D,2*pi) !== SwitchingPeriodicArray(D,4*pi) && 
      !(SwitchingPeriodicArray(D,2*pi) ≈ SwitchingPeriodicArray(D,4*pi))


@test Ad[1:2,1:1].M == Ad.M[1:2,1:1,:] && lastindex(Ad,1) == n && lastindex(Ad,2) == n

@test [[Ad Xd]; [Xd Ad]] == [[Ad;Xd] [Xd;Ad]]

@test blockdiag(Ad,Xd)[10] ≈ bldiag(Ad[10],Xd[10])   


end # pmops

end