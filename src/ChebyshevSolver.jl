
using QuadGK
include("LKvector.jl")
using .LK

struct ChebyshevSolver{isLK}
    T::Float64 #temperature
    nmax::Int64
    Tn::Vector{Float64}
    aa::Float64
    bb::Float64
    ωc::Float64
    LKeps::Float64
end

function ChebyshevSolver(T, aa, bb; kargs...)
    if :nmax in keys(kargs)
        nmax = values(kargs).nmax
    else
        nmax = 100
    end
    if :ωc in keys(kargs)
        ωc = values(kargs).nmax
    else
        ωc = 8.0
    end
    if :isLK in keys(kargs)
        isLK = values(kargs).isLK
    else
        isLK = false
    end

    if :LKeps in keys(kargs)
        LKeps = values(kargs).LKeps
    else
        LKeps = 1e-6
    end
    Tn = calc_preintegration(T, ωc, aa, bb, nmax)
    return ChebyshevSolver{isLK}(T, nmax, Tn, aa, bb, ωc, LKeps)
end

function func2(x, T, aa, bb, n)
    actb = aa * cos(x) + bb
    if actb / T > 500
        fermi = 0.0
    elseif actb / T < -500
        fermi = 1.0
    else
        fermi = 1.0 / (1 + exp(actb / T))
    end
    return -cos(n * x) * fermi
end

function calc_preintegration(T, ωc, aa, bb, nmax)
    b = (ωc - bb) / aa
    acb = acos(b)
    a = -(ωc + bb) / aa
    aca = acos(a)
    ωb = acos(-(ωc + bb) / aa)
    vec_integ = zeros(nmax + 1)

    for n = 0:nmax
        if T != 0.0
            v, err = quadgk(x -> func2(x, T, aa, bb, n), aca, acb, rtol=1e-8)
            if n == 0
                v /= 2
            end
        else
            if n == 0
                v = (ωb - ba) / 2
            else
                v = (sin(n * ωb) - sin(n * ba)) / n
            end
        end
        vec_integ[n+1] = v
    end
    return vec_integ
end

function calc_polynomials(A, x::T, left_i, right_j; nc=200) where {T}
    tempvector = Vector{T}(undef, 3)
    for i = 1:3
        tempvector[i] = zero(x)
    end
    vec_jnmm = tempvector[1]
    vec_jnm = tempvector[2]
    vec_jn = tempvector[3]
    vec_jn[right_j] = 1
    vec_ai = zeros(eltype(A), nc)

    @inbounds for nn = 0:nc-1
        if nn == 0
            vec_jnmm[right_j] = 1
        elseif nn == 1
            mul!(vec_jnmm, A, vec_jnm)
        else
            #mul!(C, A, B, α, β) -> C , A B α + C β
            mul!(vec_jnmm, A, vec_jnm, 2, -1)
        end

        vec_ai[nn+1] = vec_jnmm[left_i]
        #println(nn, "\t", vec_ai[nn+1], "\t", vec_jnmm.nonzeros)
        vec_jnm, vec_jnmm = vec_jnmm, vec_jnm
    end
    return vec_ai
end

function solve(cheby::ChebyshevSolver{false}, A::AbstractMatrix{T}, i, j) where {T}
    _, N = size(A)
    x = zeros(T, N)
    vec_a = calc_polynomials(A, x, i, j; nc=cheby.nmax)
    aa = cheby.aa
    bb = cheby.bb

    ba = acos(-bb / aa)
    omeb = acos(-(cheby.ωc + bb) / aa)

    Gij = zero(vec_a[1])
    #println(vec_a)
    #println(cheby.Tn)
    for n = 1:cheby.nmax-1
        Gij += vec_a[n+1] * cheby.Tn[n+1]
    end
    Gij += vec_a[1] * (omeb - ba) / 2
    Gij *= 2 / pi
    return Gij

end

function solve(cheby::ChebyshevSolver{false}, A::Hamiltonian{T,N,isSC,num_internal_degree,num_sites}, i, j) where {T,N,isSC,num_internal_degree,num_sites}
    _, Nt = size(A)
    x = zeros(T, Nt)
    vec_a = calc_polynomials(A.matrix, x, i, j; nc=cheby.nmax)
    aa = cheby.aa
    bb = cheby.bb

    ba = acos(-bb / aa)
    omeb = acos(-(cheby.ωc + bb) / aa)

    Gij = zero(vec_a[1])
    #println(vec_a)
    #println(cheby.Tn)
    for n = 1:cheby.nmax-1
        Gij += vec_a[n+1] * cheby.Tn[n+1]
    end
    Gij += vec_a[1] * (omeb - ba) / 2
    Gij *= 2 / pi
    return Gij

end

function solve(cheby::ChebyshevSolver{true}, A::Hamiltonian{T,N,isSC,num_internal_degree,num_sites}, i, j) where {T,N,isSC,num_internal_degree,num_sites}
    _, Nt = size(A)
    x = LKvector(T, Nt, cheby.LKeps)
    #x = zeros(T, N)
    vec_a = calc_polynomials(A.matrix, x, i, j; nc=cheby.nmax)
    aa = cheby.aa
    bb = cheby.bb

    ba = acos(-bb / aa)
    omeb = acos(-(cheby.ωc + bb) / aa)

    Gij = zero(vec_a[1])
    #println(vec_a)
    #println(cheby.Tn)
    for n = 1:cheby.nmax-1
        Gij += vec_a[n+1] * cheby.Tn[n+1]
    end
    Gij += vec_a[1] * (omeb - ba) / 2
    Gij *= 2 / pi
    return Gij

end
