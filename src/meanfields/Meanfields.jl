

struct Meanfields_solver{M,H}
    method::M
    hamiltonian::H
end

function Meanfields_solver(ham::Hamiltonian{T1,isSC}, T, ; method="RSCG", kargs...) where {T1,isSC}
    if method == "RSCG"
        println("The solver is the RSCG")
        rscg = RSCGSolver(T; kargs...)
        return Meanfields_solver{typeof(rscg),typeof(ham)}(rscg, ham)
    elseif method == "Chebyshev"
        println("The solver is the Chebyshev polynomial method")
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
        cheby = ChebyshevSolver(T, aa, bb; kargs...)
        #hamt = deepcopy(ham)
        for i = 1:N
            ham.matrix[i, i] -= bb
        end
        ham.matrix ./= aa

        return Meanfields_solver{typeof(cheby),typeof(ham)}(cheby, ham)
    else
        error("method $method is not supported yet")
    end
end

function get_hamiltonian(m::Meanfields_solver{M,Hamiltonian{T,isSC}}) where {M<:RSCGSolver,T,isSC,
    }
    return m.hamiltonian
end

function get_hamiltonian(m::Meanfields_solver{M,Hamiltonian{T,isSC}}) where {M<:ChebyshevSolver,T,isSC,
    }
    ham = deepcopy(m.hamiltonian)
    ham.matrix .*= m.method.aa
    for i = 1:N
        ham.matrix[i, i] -= m.method.bb
    end
    return ham
end
export get_hamiltonian

function calc_meanfields(m::Meanfields_solver{M,Hamiltonian{T,isSC}}, c1::FermionOP, c2::FermionOP) where {M,T,isSC,
    }
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

function calc_meanfields(m::Meanfields_solver{M,Hamiltonian{T,isSC}}, i, j) where {M,T,isSC,
    }
    Gij0 = solve(m.method, m.hamiltonian, i, j)
    return Gij0
end



function calc_greenfunction(ham::Hamiltonian{T,isSC}, zs::Vector{T2}, c1::FermionOP, c2::FermionOP) where {T,isSC,
    T2<:Number}
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

function calc_greenfunction(ham::Hamiltonian{T,isSC}, z::Number, c1::FermionOP, c2::FermionOP) where {T,isSC,
    }
    Gij = calc_greenfunction(ham, [z], c1, c2)[1]
    return Gij
end


function update_hamiltonian!(m::Meanfields_solver{M,Hamiltonian{T,isSC}}, value, c1::FermionOP, c2::FermionOP) where {M<:RSCGSolver,T,isSC,
    }

    ii = (c1.site - 1) * num_internal_degree + c1.internal_index
    jj = (c2.site - 1) * num_internal_degree + c2.internal_index
    if isSC
        ii += ifelse(c1.is_annihilation_operator, N, 0)
        jj += ifelse(c2.is_annihilation_operator, 0, N)
    else
        @assert !c1.is_annihilation_operator && c2.is_annihilation_operator "This is not C^+ C form $c1 $c2"
        #@assert !c1.is_annihilation_operator && !c2.is_annihilation_operator "This is not C^+ C^+ form $c1 $c2"
    end
    m.hamiltonian.matrix[ii, jj] = value
end

function update_hamiltonian!(m::Meanfields_solver{M,Hamiltonian{T,isSC}}, value, c1::FermionOP, c2::FermionOP) where {M<:ChebyshevSolver,T,isSC}

    ii = (c1.site - 1) * num_internal_degree + c1.internal_index
    jj = (c2.site - 1) * num_internal_degree + c2.internal_index
    if isSC
        ii += ifelse(c1.is_annihilation_operator, N, 0)
        jj += ifelse(c2.is_annihilation_operator, 0, N)
    else
        #@assert !c1.is_annihilation_operator && !c2.is_annihilation_operator "This is not C^+ C^+ form $c1 $c2"
        @assert !c1.is_annihilation_operator && c2.is_annihilation_operator "This is not C^+ C form $c1 $c2"

    end
    m.hamiltonian.matrix[ii, jj] = value / m.method.aa
    if ii == jj
        m.hamiltonian.matrix[ii, jj] -= m.method.bb / m.method.aa
    end
end
export update_hamiltonian!



export calc_meanfields, calc_greenfunction
