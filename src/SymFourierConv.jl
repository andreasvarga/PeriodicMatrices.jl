
# conversions to continuous-time Fourier function matrix
function Base.convert(::Type{FourierFunctionMatrix}, A::PeriodicSymbolicMatrix) 
    tA = convert(PeriodicFunctionMatrix,A)
    return FourierFunctionMatrix{:c,eltype(tA),Fun}(Fun(x -> tA.f(x), Fourier(0..tA.period/tA.nperiod)), Float64(tA.period), tA.nperiod)
end  

Base.convert(::Type{PeriodicSymbolicMatrix}, A::FourierFunctionMatrix) =
   convert(PeriodicSymbolicMatrix,convert(HarmonicArray,A))
