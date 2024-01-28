

struct Meanfields_solver{M,H}
    method::M
    hamiltonian::H
end

function Meanfields_solver(ham::Hamiltonian, T, ; method="RSCG", kargs...)
    if method == "RSCG"
        rscg = RSCGSolver(T; kargs...)
        return Meanfields_solver{typeof(rscg),typeof(ham)}(rscg, ham)
    else
        error("method $method is not supported yet")
    end
end

function calc_meanfields(m::Meanfields_solver{M,Hamiltonian{T,N,isSC,num_internal_degree,num_sites}}, c1::FermionOP, c2::FermionOP) where {M<:RSCGSolver,T,N,isSC,
    num_internal_degree,num_sites}
    ham = m.hamiltonian
    ii = (c1.site - 1) * num_internal_degree + c1.internal_index
    jj = (c2.site - 1) * num_internal_degree + c2.internal_index
    if isSC
        ii += ifelse(c1.is_annihilation_operator, N, 0)
        jj += ifelse(c2.is_annihilation_operator, 0, N)
    else
        @assert !c1.is_annihilation_operator && !c2.is_annihilation_operator "This is not C^+ C^+ form $c1 $c2"
    end

    Gij0 = solve(m.method, m.hamiltonian, jj, ii)
    return Gij0
end

function calc_meanfields(m::Meanfields_solver{M,Hamiltonian{T,N,isSC,num_internal_degree,num_sites}}, i, j) where {M<:RSCGSolver,T,N,isSC,
    num_internal_degree,num_sites}
    Gij0 = solve(m.method, m.hamiltonian, i, j)
    return Gij0
end

function calc_greenfunction(ham::Hamiltonian{T,N,isSC,num_internal_degree,num_sites}, zs::Vector{T2}, c1::FermionOP, c2::FermionOP) where {T,N,isSC,
    num_internal_degree,num_sites,T2<:Number}
    ii = (c1.site - 1) * num_internal_degree + c1.internal_index
    jj = (c2.site - 1) * num_internal_degree + c2.internal_index
    if isSC
        ii += ifelse(c1.is_annihilation_operator, N, 0)
        jj += ifelse(c2.is_annihilation_operator, 0, N)
    else
        @assert !c1.is_annihilation_operator && !c2.is_annihilation_operator "This is not C^+ C^+ form $c1 $c2"
    end
    Gijs = greensfunctions(jj, ii, zs, ham)
    return Gijs
end

function calc_greenfunction(ham::Hamiltonian{T,N,isSC,num_internal_degree,num_sites}, z::Number, c1::FermionOP, c2::FermionOP) where {T,N,isSC,
    num_internal_degree,num_sites}
    Gij = calc_greenfunction(ham, [z], c1, c2)[1]
    return Gij
end




export calc_meanfields, calc_greenfunction