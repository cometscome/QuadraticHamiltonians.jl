using RSCG
struct ContourProjector
    α::Float64
    ρ::Float64 #radious
    γ::Float64 #center
    Nq::Int64
    vec_w::Vector{ComplexF64} #weights
    vec_zP::Vector{ComplexF64} #frequencies
    vec_zQ::Vector{ComplexF64} #frequencies
    eps::Float64

    function ContourProjector(Emin, α, Nq,eps)
        ρ = abs(Emin / 2)
        γ = Emin / 2

        vec_zP = zeros(ComplexF64, Nq)
        vec_zQ = zeros(ComplexF64, Nq)
        vec_w = zeros(ComplexF64, Nq)

        for j in 1:Nq
            θj = (2π / Nq) * (j - 1 / 2)
            vec_zP[j] = γ + ρ * (cos(θj) + im * α * sin(θj))
            vec_zQ[j] = -γ + ρ * (cos(θj) + im * α * sin(θj))
            vec_w[j] = α * cos(θj) + im * sin(θj)
        end
        return new(α, ρ, γ, Nq, vec_w, vec_zP, vec_zQ,eps)
    end
end

function apply_P!(projector::Projector{M,H}, Px, x, ϵ=0.0) where {M<:ContourProjector,H}
    method = projector.method!
    vec_x = greensfunctions_col(method.vec_zP, projector.hamiltonian.matrix, x, eps=method.eps)

    dim = get_dim(projector.hamiltonian)
    Px .= 0
    for i = 1:dim
        for j in 1:method.Nq
            Px[i] += method.ρ * method.vec_w[j] * vec_x[j][i] / method.Nq
        end
    end
end

function apply_Q!(projector::Projector{M,H}, Px, x, ϵ=0.0) where {M<:ContourProjector,H}
    method = projector.method!
    vec_x = greensfunctions_col(method.vec_zQ, projector.hamiltonian.matrix, x, eps=method.eps)

    dim = get_dim(projector.hamiltonian)
    Px .= 0
    for i = 1:dim
        for j in 1:method.Nq
            Px[i] += method.ρ * method.vec_w[j] * vec_x[j][i] / method.Nq
        end
    end
end

