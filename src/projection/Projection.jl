using LinearAlgebra
struct Projector{M,H}
    method!::M
    hamiltonian::H
end
export Projector

function Projector(ham::Hamiltonian{T1,N,isSC,num_internal_degree,num_sites}, ; method="Exact", kargs...) where {T1,N,isSC,num_internal_degree,num_sites}
    if method == "Exact"
        println("The exact diagonalization method is used")
        exact = ExactProjector(ham)
        return Projector{typeof(exact),typeof(ham)}(exact, ham)
    elseif method == "KPM"
        dim = get_dim(ham)
        println("The solver is the Kernel Polynomial Method")
        if :aa in keys(kargs)
            aa = values(kargs).aa
        else
            aa = 10.0
        end
        if :bb in keys(kargs)
            bb = values(kargs).bb
        else
            bb = 0.0
        end
        kpm = KPMProjector(T, aa, bb; kargs...)
        #hamt = deepcopy(ham)
        for i = 1:dim
            ham.matrix[i, i] -= bb
        end
        ham.matrix ./= aa

        return Projector{typeof(kpm),typeof(ham)}(kpm, ham)
    elseif method == "Contour"
        dim = get_dim(ham)
        println("The solver is the contour integral Method")
        if :Emin in keys(kargs)
            Emin = values(kargs).Emin
        else
            Emin = -10
        end





        if :α in keys(kargs)
            α = values(kargs).α
        else
            α = 0.1
        end

        if :Nq in keys(kargs)
            Nq = values(kargs).Nq
        else
            Nq = 100
        end





        p = ContourProjector(Emin, α, Nq)
        return Projector{typeof(p),typeof(ham)}(p, ham)
    else
        error("method $method is not supported yet")
    end
end

struct ExactProjector{T1,T2}
    eigenvalues::Vector{T1}
    eigenvectors::Matrix{T2}

    function ExactProjector(ham::Hamiltonian{T1,N,isSC,num_internal_degree,num_sites}) where {T1,N,isSC,num_internal_degree,num_sites}
        evals, evecs = eigen(Matrix(ham.matrix))
        return new{eltype(evals),eltype(evecs)}(evals, evecs)
    end
end
export ExactProjector

function (m::ExactProjector)(x; ϵ=0.0)
    Px = zero(x)
    m(Px, x, ϵ)
    return Px
end

function (m::ExactProjector)(Px, x; ϵ=0.0)
    Px .= 0
    for (n, energy) in enumerate(m.eigenvalues)
        if energy < ϵ
            a = dot(m.eigenvectors[:, n], x)
            Px .+= a * m.eigenvectors[:, n]
        end
    end
end

function apply_P!(projector::Projector{M,H}, Px, x, ϵ=0.0) where {M<:ExactProjector,H}
    projector.method!(Px, x; ϵ)
end


function apply_Q!(projector::Projector{M,H}, Px, x, ϵ=0.0) where {M<:ExactProjector,H}
    projector.method!(Px, x; ϵ)
    for (i, xi) in enumerate(x)
        Px[i] = xi - Px[i]
    end
end


function make_PQ(m::ExactProjector)
    v = m.eigenvectors
    N, _ = size(v)
    e = m.eigenvalues
    P = zero(v)
    Q = zero(v)
    N = length(e)
    for i = 1:N
        vv = v[:, i]
        if e[i] < 0
            P += vv * vv'
        else
            Q += vv * vv'
        end
    end
    return P, Q
end



function make_C(p::Projector{M,H}, isite, Nx; debugmode=false) where {M,T,N,isSC,num_internal_degree,num_sites,H<:Hamiltonian{T,N,isSC,num_internal_degree,num_sites}}
    X = make_X(p.hamiltonian, Nx)
    Y = make_Y(p.hamiltonian, Nx)
    n = get_dim(p.hamiltonian)
    e = zeros(T, n)


    if debugmode
        C = zero(T)
        P, Q = make_PQ(p.method!)
        xQ = Q * X * P
        yP = P * Y * Q
        Cmatrix = xQ * yP
        for ispin = 1:num_internal_degree
            i = (isite - 1) * num_internal_degree + ispin
            Ci = -4pi * imag(Cmatrix[i, i]) - 4pi * imag(Cmatrix[i+N, i+N])
            C += Ci
        end
        println(C)
    end
    C = zero(T)

    for ispin = 1:num_internal_degree
        i = (isite - 1) * num_internal_degree + ispin
        e .= 0
        e[i] = 1
        tempvector1 = zero(e)
        tempvector2 = zero(e)

        apply_Q!(p, tempvector1, e)
        mul!(tempvector2, Y, tempvector1)
        apply_P!(p, tempvector1, tempvector2)
        mul!(tempvector2, X, tempvector1)
        apply_Q!(p, tempvector1, tempvector2)

        Ci = -4pi * imag(tempvector1[i])
        C += Ci
        #println(Ci)

        e .= 0
        e[i+num_internal_degree*num_sites] = 1
        apply_Q!(p, tempvector1, e)
        mul!(tempvector2, Y, tempvector1)
        apply_P!(p, tempvector1, tempvector2)
        mul!(tempvector2, X, tempvector1)
        apply_Q!(p, tempvector1, tempvector2)
        Ci = -4pi * imag(tempvector1[i+num_internal_degree*num_sites])
        C += Ci

        #println(Ci)

    end
    #println(C)
    #error("d")

    return C


end
export make_C
