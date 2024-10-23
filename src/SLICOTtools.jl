module SLICOTtools
# Collection of wrappers extracted from SLICOTMath.jl (generated by Ralph Smith)
# based on the SLICOT_jll library (created by the courtesy of Ralph Smith)
using SLICOT_jll
using LinearAlgebra
using LinearAlgebra: BlasInt

export mb03vd!, mb03vy!, mb03bd!, mb03wd!, mb03kd!

function chkargsok(ret::BlasInt)
    if ret < 0
        throw(ArgumentError("invalid argument #$(-ret) to SLICOT call"))
    end
end

const BlasBool = BlasInt

"""
    mb03vd!(n::Integer, p::Integer, ilo::Integer, ihi::Integer, A::Array{Float64, 3}, tau::AbstractMatrix{Float64}) -> info::Int64

Reduce a product of `p` real general matrices `A = A_1*A_2*...*A_p`
to upper Hessenberg form, `H = H_1*H_2*...*H_p`, where `H_1` is
upper Hessenberg, and `H_2`, ..., `H_p` are upper triangular, by using
orthogonal similarity transformations on `A`,

        Q_1' * A_1 * Q_2 = H_1,
        Q_2' * A_2 * Q_3 = H_2,
               ...
        Q_p' * A_p * Q_1 = H_p.

The matrices `A_1`, `A_2`, ..., `A_p` are contained in the 3-dimensional array `A`. 
The resulting `H_1`, `H_2`, ..., `H_p` and `Q_1`, `Q_2`, ..., `Q_p` overwrite 
`A_1`, `A_2`, ..., `A_p` in `A` and the array `tau`.   

See the SLICOT documentation of `MB03VD` for details.
"""
function mb03vd!(n::Integer, p::Integer, ilo::Integer, ihi::Integer,
    a::Array{Float64,3}, tau::AbstractMatrix{Float64})

    lda1 = max(1,stride(a,2))
    lda2 = max(1,stride(a,3)÷lda1)
    ldtau = max(1,stride(tau,2))
    info = Ref{BlasInt}()
    dwork = Vector{Float64}(undef, n)

    ccall((:mb03vd_, libslicot), Cvoid, (Ref{BlasInt}, Ref{BlasInt},
            Ref{BlasInt}, Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt},
            Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt}, Ptr{Float64},
            Ptr{BlasInt}), n, p, ilo, ihi, a, lda1, lda2, tau,
            ldtau, dwork, info)
    chkargsok(info[])

    return info[]
end

"""
     mb03vy!(n::Integer, p::Integer, ilo::Integer, ihi::Integer, A::Array{Float64, 3}, tau::AbstractMatrix{Float64}) -> info::Int64

Generate the real orthogonal matrices `Q_1`, `Q_2`, ..., `Q_p`,
which are defined as the product of `ihi-ilo` elementary reflectors
of order `n`, as returned in `A_1`, `A_2`, ..., `A_p` by `mb03vd!`:

     Q_j = H_j(ilo) H_j(ilo+1) . . . H_j(ihi-1).

The 3-dimensional arrays `A` and `tau` contains the information on the employed
elementary reflectors. The resulting `Q_1`, `Q_2`, ..., `Q_p` overwrite `A_1`, `A_2`, ..., `A_p`. 

See the SLICOT documentation of `MB03VY` for details.

"""
function mb03vy!(n::Integer, p::Integer, ilo::Integer, ihi::Integer,
    a::Array{Float64,3}, tau::AbstractMatrix{Float64})

    lda1 = max(1,stride(a,2))
    lda2 = max(1,stride(a,3)÷lda1)
    ldtau = max(1,stride(tau,2))
    info = Ref{BlasInt}()
    ldwork = BlasInt(-1)
    dwork = Vector{Float64}(undef, 1)

    local jlres
    for iwq in 1:2
        ccall((:mb03vy_, libslicot), Cvoid, (Ref{BlasInt}, Ref{BlasInt},
            Ref{BlasInt}, Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt},
            Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt}, Ptr{Float64},
            Ref{BlasInt}, Ptr{BlasInt}), n, p, ilo, ihi, a, lda1,
            lda2, tau, ldtau, dwork, ldwork, info)
        chkargsok(info[])
        if iwq == 1
            ldwork = BlasInt(real(dwork[1]))
            resize!(dwork, ldwork)
        end
    end

    return info[]
end


"""
    mb03vw!(compq::AbstractChar, triu::AbstractChar, qind::AbstractVector{Int64}, k::Integer, n::Integer, h::Integer, 
            ilo::Integer, ihi::Integer, s::AbstractVector{Int64}, a::Array{Float64, 3}, q::Array{Float64, 3}, 
            liwork::Integer, ldwork::Integer) -> info::Int64
            
    mb03vw!(compq::AbstractChar, triu::AbstractChar, qind::AbstractVector{Int64}, k::Integer, n::Integer, h::Integer, 
            ilo::Integer, ihi::Integer, s::AbstractVector{Int64}, a::Array{Float64, 3}, q::Array{Float64, 3}, 
            iwork::AbstractVector{Int64}, dwork::AbstractVector{Float64}) -> info::Int64

Reduce the generalized matrix product

              s[1]           s[2]                 s[k]
      A[:,:,1]     * A[:,:,2]     * ... * A[:,:,k]

to upper Hessenberg-triangular form, where A is N-by-N-by-K and S
is the signature array with values 1 or -1. The H-th matrix of A
is reduced to upper Hessenberg form while the other matrices are
triangularized. 

If `compq = 'U'` or `compq = 'I'`, then the orthogonal factors are
computed and stored in the array `Q` so that for `s[i] = 1`,

                    T
        Q[:,:,i](in)   A[:,:,i](in)   Q[:,:,mod(i,k)+1](in)
                                                            T 
    =   Q[:,:,i](out)  A[:,:,i](out)  Q[:,:,mod(i,k)+1](out),

and for `s[i] = -1`,

                             T
        Q[:,:,mod(i,k)+1](in)   A[:,:,i](in)   Q[:,:,i](in)
                                                            T 
    =   Q[:,:,mod(i,k)+1](out)  A[:,:,i](out)  Q[:,:,i](out).

A partial generation of the orthogonal factors can be realized
via the array `qind`.

If `triu = 'N'` only matrices with negative signature are reduced to upper
triangular form in the first stage of the algorithm. 
If `triu = 'A'` all possible `n-1` matrices with negative signature are reduced. 

See the SLICOT documentation of `MB03VW` for details.
"""
function mb03vw!(compq::AbstractChar, qind::AbstractVector{BlasInt}, atriu::AbstractChar, 
    n::Integer, k::Integer, h::Integer, ilo::Integer, ihi::Integer,
    s::AbstractVector{BlasInt}, a::Array{Float64,3},
    q::Array{Float64,3}, liwork::Integer, ldwork::Integer)

    lda1 = max(1,stride(a,2))
    lda2 = max(1,stride(a,3)÷lda1)
    ldq1 = max(1,stride(q,2))
    ldq2 = max(1,stride(q,3)÷ldq1)
    info = Ref{BlasInt}()
    iwork = Vector{BlasInt}(undef, liwork)
    dwork = Vector{Float64}(undef, ldwork)
    hi = (h < 0 || h > k) ? 1 : h

    ccall((:mb03vw_, libslicot), Cvoid, (Ref{UInt8}, Ptr{BlasInt}, 
            Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, 
            Ref{BlasInt}, Ref{BlasInt}, Ptr{BlasInt}, Ptr{Float64}, 
            Ref{BlasInt}, Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt}, 
            Ref{BlasInt}, Ptr{BlasInt}, Ref{BlasInt},
            Ptr{Float64}, Ref{BlasInt}, Ptr{BlasInt}, Clong, Clong), 
            compq, qind, atriu, n, k, hi,
            ilo, ihi, s, a, lda1, lda2, q, ldq1, ldq2, 
            iwork, liwork, dwork, ldwork, info, 1, 1)
    chkargsok(info[])

    return info[]
end
function mb03vw!(compq::AbstractChar, qind::AbstractVector{BlasInt}, atriu::AbstractChar, 
    n::Integer, k::Integer, h::Integer, ilo::Integer, ihi::Integer,
    s::AbstractVector{BlasInt}, a::Array{Float64,3},
    q::Array{Float64,3}, iwork::AbstractVector{BlasInt}, dwork::AbstractVector{Float64})

    lda1 = max(1,stride(a,2))
    lda2 = max(1,stride(a,3)÷lda1)
    ldq1 = max(1,stride(q,2))
    ldq2 = max(1,stride(q,3)÷ldq1)
    info = Ref{BlasInt}()
    liwork = length(iwork)
    ldwork = length(dwork)
    hi = (h < 0 || h > k) ? 1 : h

    ccall((:mb03vw_, libslicot), Cvoid, (Ref{UInt8}, Ptr{BlasInt}, 
            Ref{UInt8}, Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, 
            Ref{BlasInt}, Ref{BlasInt}, Ptr{BlasInt}, Ptr{Float64}, 
            Ref{BlasInt}, Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt}, 
            Ref{BlasInt}, Ptr{BlasInt}, Ref{BlasInt},
            Ptr{Float64}, Ref{BlasInt}, Ptr{BlasInt}, Clong, Clong), 
            compq, qind, atriu, n, k, hi,
            ilo, ihi, s, a, lda1, lda2, q, ldq1, ldq2, 
            iwork, liwork, dwork, ldwork, info, 1, 1)
    chkargsok(info[])

    return info[]
end

"""
    mb03bd!(job::AbstractChar, defl::AbstractChar, compq::AbstractChar, qind::AbstractVector{Int64}, k::Integer, n::Integer, h::Integer, 
            ilo::Integer, ihi::Integer, s::AbstractVector{Int64}, a::Array{Float64, 3}, q::Array{Float64, 3}, alphar::AbstractVector{Float64}, 
            alphai::AbstractVector{Float64}, beta::AbstractVector{Float64}, scal::AbstractVector{Int64}, 
            liwork::Integer, ldwork::Integer) -> (info::Int64, iwarn::Int64)

    mb03bd!(job::AbstractChar, defl::AbstractChar, compq::AbstractChar, qind::AbstractVector{Int64}, k::Integer, n::Integer, h::Integer, 
            ilo::Integer, ihi::Integer, s::AbstractVector{Int64}, a::Array{Float64, 3}, q::Array{Float64, 3}, alphar::AbstractVector{Float64}, 
            alphai::AbstractVector{Float64}, beta::AbstractVector{Float64}, scal::AbstractVector{Int64}, 
            iwork::AbstractVector{Int64}, dwork::AbstractVector{Float64}) -> (info::Int64, iwarn::Int64)

Find the eigenvalues of the generalized matrix product

              s[1]           s[2]                 s[k]
      A[:,:,1]     * A[:,:,2]     * ... * A[:,:,k]

where `A[:,:,h]` is upper Hessenberg and `A[:,:,i]`, `i <> h`, is upper
triangular, using a double-shift version of the periodic
QZ method. In addition, `A` may be reduced to periodic Schur form:
`A[:,:,h]` is upper quasi-triangular and all the other factors
`A[:,:,i]` are upper triangular. Optionally, the 2-by-2 triangular
matrices corresponding to 2-by-2 diagonal blocks in `A[:,:,h]`
are so reduced that their product is a 2-by-2 diagonal matrix.

If `compq = 'U'` or `compq = 'I'`, then the orthogonal factors are
computed and stored in the array `Q` so that for `s[i] = 1`,

                    T
        Q[:,:,i](in)   A[:,:,i](in)   Q[:,:,mod(i,k)+1](in)
                                                            T 
    =   Q[:,:,i](out)  A[:,:,i](out)  Q[:,:,mod(i,k)+1](out),

and for `s[i] = -1`,

                             T
        Q[:,:,mod(i,k)+1](in)   A[:,:,i](in)   Q[:,:,i](in)
                                                            T 
    =   Q[:,:,mod(i,k)+1](out)  A[:,:,i](out)  Q[:,:,i](out).

A partial generation of the orthogonal factors can be realized
via the array `qind`.

See the SLICOT documentation of `MB03BD` for details.
"""
function mb03bd!(job::AbstractChar, defl::AbstractChar,
    compq::AbstractChar, qind::AbstractVector{BlasInt}, k::Integer,
    n::Integer, h::Integer, ilo::Integer, ihi::Integer,
    s::AbstractVector{BlasInt}, a::Array{Float64,3},
    q::Array{Float64,3},
    alphar::AbstractVector{Float64}, alphai::AbstractVector{Float64},
    beta::AbstractVector{Float64}, scal::AbstractVector{BlasInt},
    liwork::Integer, ldwork::Integer)

    lda1 = max(1,stride(a,2))
    lda2 = max(1,stride(a,3)÷lda1)
    ldq1 = max(1,stride(q,2))
    ldq2 = max(1,stride(q,3)÷ldq1)
    info = Ref{BlasInt}()
    iwarn = Ref{BlasInt}()
    iwork = Vector{BlasInt}(undef, liwork)
    dwork = Vector{Float64}(undef, ldwork)
    ccall((:mb03bd_, libslicot), Cvoid, (Ref{UInt8}, Ref{UInt8},
            Ref{UInt8}, Ptr{BlasInt}, Ref{BlasInt}, Ref{BlasInt},
            Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ptr{BlasInt},
            Ptr{Float64}, Ref{BlasInt}, Ref{BlasInt}, Ptr{Float64},
            Ref{BlasInt}, Ref{BlasInt}, Ptr{Float64}, Ptr{Float64},
            Ptr{Float64}, Ptr{BlasInt}, Ptr{BlasInt}, Ref{BlasInt},
            Ptr{Float64}, Ref{BlasInt}, Ptr{BlasInt}, Ptr{BlasInt},
            Clong, Clong, Clong), job, defl, compq, qind, k, n, h,
            ilo, ihi, s, a, lda1, lda2, q, ldq1, ldq2, alphar,
            alphai, beta, scal, iwork, liwork, dwork, ldwork, iwarn,
            info, 1, 1, 1)
    chkargsok(info[])

    return info[], iwarn[]
end
function mb03bd!(job::AbstractChar, defl::AbstractChar,
    compq::AbstractChar, qind::AbstractVector{BlasInt}, k::Integer,
    n::Integer, h::Integer, ilo::Integer, ihi::Integer,
    s::AbstractVector{BlasInt}, a::Array{Float64,3},
    q::Array{Float64,3},
    alphar::AbstractVector{Float64}, alphai::AbstractVector{Float64},
    beta::AbstractVector{Float64}, scal::AbstractVector{BlasInt},
    iwork::AbstractVector{BlasInt}, dwork::AbstractVector{Float64})

    lda1 = max(1,stride(a,2))
    lda2 = max(1,stride(a,3)÷lda1)
    ldq1 = max(1,stride(q,2))
    ldq2 = max(1,stride(q,3)÷ldq1)
    info = Ref{BlasInt}()
    iwarn = Ref{BlasInt}()
    liwork = length(iwork)
    ldwork = length(dwork)
    ccall((:mb03bd_, libslicot), Cvoid, (Ref{UInt8}, Ref{UInt8},
            Ref{UInt8}, Ptr{BlasInt}, Ref{BlasInt}, Ref{BlasInt},
            Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ptr{BlasInt},
            Ptr{Float64}, Ref{BlasInt}, Ref{BlasInt}, Ptr{Float64},
            Ref{BlasInt}, Ref{BlasInt}, Ptr{Float64}, Ptr{Float64},
            Ptr{Float64}, Ptr{BlasInt}, Ptr{BlasInt}, Ref{BlasInt},
            Ptr{Float64}, Ref{BlasInt}, Ptr{BlasInt}, Ptr{BlasInt},
            Clong, Clong, Clong), job, defl, compq, qind, k, n, h,
            ilo, ihi, s, a, lda1, lda2, q, ldq1, ldq2, alphar,
            alphai, beta, scal, iwork, liwork, dwork, ldwork, iwarn,
            info, 1, 1, 1)
    chkargsok(info[])

    return info[], iwarn[]
end


function mb03bz!(job::AbstractChar, compq::AbstractChar, k::Integer,
    n::Integer, ilo::Integer, ihi::Integer,
    s::AbstractVector{BlasInt}, a::Array{ComplexF64,3}, q::Array{ComplexF64,3},
    alpha::AbstractVector{ComplexF64},
    beta::AbstractVector{ComplexF64}, scal::AbstractVector{BlasInt},
    ldwork::Integer, lzwork::Integer)
"""
    mb03bz!(job::AbstractChar, compq::AbstractChar, k::Integer, n::Integer, 
            ilo::Integer, ihi::Integer, s::AbstractVector{Int64}, a::Array{ComplexF64, 3}, 
            q::Array{ComplexF64, 3}, alpha::AbstractVector{ComplexF64}, beta::AbstractVector{ComplexF64}, 
            scal::AbstractVector{Int64}, ldwork::Integer, lzwork::Integer) -> info::Int64

"""

    lda1 = max(1,stride(a,2))
    lda2 = max(1,stride(a,3)÷lda1)
    ldq1 = max(1,stride(q,2))
    ldq2 = max(1,stride(q,3)÷ldq1)
    info = Ref{BlasInt}()
    dwork = Vector{Float64}(undef, ldwork)
    zwork = Vector{ComplexF64}(undef, lzwork)

    ccall((:mb03bz_, libslicot), Cvoid, (Ref{UInt8}, Ref{UInt8},
            Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt},
            Ptr{BlasInt}, Ptr{ComplexF64}, Ref{BlasInt},
            Ref{BlasInt}, Ptr{ComplexF64}, Ref{BlasInt},
            Ref{BlasInt}, Ptr{ComplexF64}, Ptr{ComplexF64},
            Ptr{BlasInt}, Ptr{Float64}, Ref{BlasInt},
            Ptr{ComplexF64}, Ref{BlasInt}, Ptr{BlasInt}, Clong,
            Clong), job, compq, k, n, ilo, ihi, s, a, lda1, lda2, q,
            ldq1, ldq2, alpha, beta, scal, dwork, ldwork, zwork,
            lzwork, info, 1, 1)
    chkargsok(info[])

    return info[]
end
"""
    mb03wd!(job::AbstractChar, compz::AbstractChar, n::Integer, p::Integer, 
            ilo::Integer, ihi::Integer, iloz::Integer, ihiz::Integer,  h::Array{Float64, 3}, z::Array{Float64, 3}, 
            wr::AbstractVector{Float64}, wi::AbstractVector{Float64}, ldwork::Integer) -> info::Int64

Compute the Schur decomposition and the eigenvalues of a
product of matrices, H = `H_1`*`H_2`*...*`H_p`, with `H_1` an upper
Hessenberg matrix and `H_2`, ..., `H_p` upper triangular matrices,
without evaluating the product. Specifically, the matrices Z_i
are computed, such that

        `Z_1' * H_1 * Z_2 = T_1,`
        `Z_2' * H_2 * Z_3 = T_2,`
               `...`
        `Z_p' * H_p * Z_1 = T_p,`

where `T_1` is in real Schur form, and `T_2`, ..., `T_p` are upper
triangular.

The routine works primarily with the Hessenberg and triangular
submatrices in rows and columns ILO to IHI, but optionally applies
the transformations to all the rows and columns of the matrices
H_i, i = 1,...,p. The transformations can be optionally
accumulated.

See the SLICOT documentation of `MB03WD` for details.
"""
function mb03wd!(job::AbstractChar, compz::AbstractChar, n::Integer,
    p::Integer, ilo::Integer, ihi::Integer, iloz::Integer,
    ihiz::Integer, h::Array{Float64,3}, z::Array{Float64,3},
    wr::AbstractVector{Float64}, wi::AbstractVector{Float64},
    ldwork::Integer)

    ldh1 = max(1,stride(h,2))
    ldh2 = max(1,stride(h,3)÷ldh1)
    ldz1 = max(1,stride(z,2))
    ldz2 = max(1,stride(z,3)÷ldz1)
    info = Ref{BlasInt}()
    dwork = Vector{Float64}(undef, ldwork)

    ccall((:mb03wd_, libslicot), Cvoid, (Ref{UInt8}, Ref{UInt8},
            Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt},
            Ref{BlasInt}, Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt},
            Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt}, Ref{BlasInt},
            Ptr{Float64}, Ptr{Float64}, Ptr{Float64}, Ref{BlasInt},
            Ptr{BlasInt}, Clong, Clong), job, compz, n, p, ilo, ihi,
            iloz, ihiz, h, ldh1, ldh2, z, ldz1, ldz2, wr, wi, dwork,
            ldwork, info, 1, 1)
    chkargsok(info[])

    return info[]
end
"""
    mb03kd!(compq::AbstractChar, strong::AbstractChar, k::Integer, nc::Integer, kschur::Integer, n::AbstractVector{Int64}, ni::AbstractVector{Int64}, 
            s::AbstractVector{Int64}, select::AbstractVector{BlasInt}, t::AbstractVector{Float64}, ldt::AbstractVector{Int64}, ixt::AbstractVector{Int64}, 
            q::AbstractVector{Float64}, ldq::AbstractVector{Int64}, ixq::AbstractVector{Int64}, tol::Float64, ldwork::Integer) -> (m::Int64, info::Int64)

Reorder the diagonal blocks of the formal matrix product

     T22_k^s[k] * T22_k-1^s[k-1] * ... * T22_1^s[1],                (1)

of length `k`, in the generalized periodic Schur form,

              [  T11_i  T12_i  T13_i  ]
        T_i = [    0    T22_i  T23_i  ],    i = 1, ..., k,          (2)
              [    0      0    T33_i  ]

where

  - the submatrices `T11_i` are `ni(i+1)-by-ni(i)`, if `s[i] = 1`, or
    `ni(i)-by-ni(i+1)`, if `s[i] = -1`, and contain dimension-induced infinite eigenvalues,

  - the submatrices `T22_i` are `nc-by-nc` and contain core eigenvalues,
    which are generically neither zero nor infinite,

  - the submatrices `T33_i` contain dimension-induced zero eigenvalues,

such that the `m` selected eigenvalues pointed to by the integer
vector `select` end up in the leading part of the matrix sequence
`T22_i`.

Given that `n[i] = n[i+1]` for all `i` where s`[i] = -1`, the `T11_i` are
void and the first `m` columns of the updated orthogonal
transformation matrix sequence `Q_1, ..., Q_k` span a periodic
deflating subspace corresponding to the same eigenvalues.


If `compq = 'U'` or `compq = 'I'`, then the orthogonal factors are
computed and stored in the array `Q` so that for `s[i] = 1`,

               T
        Q_i(in)  T_i(in) Q_(mod(i,k)+1)(in)
                                              T 
    =   Q_i(out) T_i(out)  Q_(mod(i,k)+1)(out),

and for `s[i] = -1`,

                          T
        Q_(mod(i,k)+1)(in)   T_i(in)   Q_i(in)
                                               T 
    =   Q_(mod(i,k)+1)(out)  T_i(out)  Q_i(out).


See the SLICOT documentation of `MB03KD` for details.
"""
function mb03kd!(compq::AbstractChar, strong::AbstractChar, k::Integer, nc::Integer, kschur::Integer, 
    n::AbstractVector{BlasInt}, ni::AbstractVector{BlasInt}, 
    s::AbstractVector{BlasInt}, select::AbstractVector{BlasInt}, 
    t::AbstractVector{Float64}, ldt::AbstractVector{BlasInt}, ixt::AbstractVector{BlasInt}, 
    q::AbstractVector{Float64}, ldq::AbstractVector{BlasInt}, ixq::AbstractVector{BlasInt}, 
    tol::Float64, ldwork::Integer)

    whichq = Vector{BlasInt}(undef, 1)
    m = Ref{BlasInt}()
    info = Ref{BlasInt}()
    iwork = Vector{BlasInt}(undef, 4*k)
    dwork = Vector{Float64}(undef, ldwork)

    ccall((:mb03kd_, libslicot), Cvoid, (Ref{UInt8}, Ptr{BlasInt}, Ref{UInt8},
            Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, 
            Ptr{BlasInt}, Ptr{BlasInt}, Ptr{BlasInt}, Ptr{BlasInt}, 
            Ptr{Float64}, Ptr{BlasInt}, Ptr{BlasInt}, Ptr{Float64}, Ptr{BlasInt}, Ptr{BlasInt}, 
            Ptr{BlasInt},  Ref{Float64}, Ptr{BlasInt}, Ptr{Float64}, Ref{BlasInt}, Ptr{BlasInt}, 
            Clong, Clong), compq, whichq, strong, k, nc, kschur, n, ni, s, select, 
            t, ldt, ixt, q, ldq, ixq, m, tol, iwork, dwork, ldwork, info, 1, 1)
    chkargsok(info[])

    return m[], info[]
end


end # module
