module QuadraticHamiltonians
export QuadraticOPs, FermionOP, Hamiltonian,
    construct_matrix
# Write your package code here.

include("RSCGSolver.jl")
using .Mod_RSCGSolver
import .Mod_RSCGSolver: RSCGSolver
export RSCGSolver


struct FermionOP{T}
    site::Int64
    internal_index::Int64
    is_annihilation_operator::Bool
    value::T
    function FermionOP(site, internal_index, is_annihilation_operator, value)
        return new{typeof(value)}(site, internal_index, is_annihilation_operator, value)
    end
end
function FermionOP(site, internal_index)
    return FermionOP(site, internal_index, true, 1)
end
function FermionOP(site)
    return FermionOP(site, 1, true, 1)
end

Base.adjoint(c::FermionOP) = FermionOP(c.site, c.internal_index, !c.is_annihilation_operator, c.value)

struct QuadraticOPs{T}
    operators::Vector{Tuple{FermionOP{T},FermionOP{T}}}
    values::Vector{T}
end

function QuadraticOPs(T::DataType=Float64)
    return QuadraticOPs{T}(Tuple{FermionOP{T},FermionOP{T}}[], T[])
end

function find_operator(h1::QuadraticOPs, op::Tuple{FermionOP{T},FermionOP{T}}) where {T}
    for (i, operator) in enumerate(h1.operators)
        c1 = operator[1]
        if c1.site == op[1].site && c1.internal_index == op[1].internal_index &&
           c1.is_annihilation_operator == op[1].is_annihilation_operator
            c2 = operator[2]
            if c2.site == op[2].site && c2.internal_index == op[2].internal_index &&
               c2.is_annihilation_operator == op[2].is_annihilation_operator
                return i
            end
        end
    end
    return 0
end

function Base.:*(a::T, c::FermionOP) where {T<:Number}
    return FermionOP(c.site, c.internal_index, c.is_annihilation_operator, a * c.value)
end

function Base.:*(c::FermionOP, a::T) where {T<:Number}
    return FermionOP(c.site, c.internal_index, c.is_annihilation_operator, a * c.value)
end

function Base.:*(c1::FermionOP, c2::FermionOP)
    value = c1.value * c2.value

    operators = [(FermionOP(c1.site, c1.internal_index, c1.is_annihilation_operator, one(value)),
        FermionOP(c2.site, c2.internal_index, c2.is_annihilation_operator, one(value)))]
    values = [value]
    return QuadraticOPs(operators, values)
end

function Base.:+(h1::QuadraticOPs{T1}, h2::QuadraticOPs{T2}) where {T1,T2}
    T = T1
    if T2 <: Complex
        T = ComplexF64
    end
    operators = Tuple{FermionOP{T},FermionOP{T}}[]
    values = T[]
    for (i, operator) in enumerate(h1.operators)
        c1 = operator[1]
        c2 = operator[2]
        value = h1.values[i]
        push!(operators, (FermionOP(c1.site, c1.internal_index, c1.is_annihilation_operator, one(T)),
            FermionOP(c2.site, c2.internal_index, c2.is_annihilation_operator, one(T))))
        push!(values, value)
    end

    for (i, operator) in enumerate(h2.operators)
        c1 = operator[1]
        c2 = operator[2]
        value = h2.values[i]
        op = (FermionOP(c1.site, c1.internal_index, c1.is_annihilation_operator, one(T)),
            FermionOP(c2.site, c2.internal_index, c2.is_annihilation_operator, one(T)))
        position = find_operator(h1, op)
        if position == 0 #not found
            push!(operators, op)
            push!(values, value)
        else
            values[position] += value
        end
    end
    return QuadraticOPs(operators, values)
end

function Base.adjoint(h1::QuadraticOPs{T1}) where {T1}
    T = T1
    operators = Tuple{FermionOP{T},FermionOP{T}}[]
    values = T[]
    for (i, operator) in enumerate(h1.operators)
        c1 = operator[1]
        c2 = operator[2]
        value = h1.values[i]'
        push!(operators, (
            FermionOP(c2.site, c2.internal_index, !c2.is_annihilation_operator, one(value)),
            FermionOP(c1.site, c1.internal_index, !c1.is_annihilation_operator, one(value))))
        push!(values, value)
    end
    return QuadraticOPs(operators, values)
end


function Base.display(h::QuadraticOPs)
    cdagcstring = ""
    for (i, operator) in enumerate(h.operators)
        c1 = operator[1]

        #println(c1)
        cname1 = "C"
        cname1 *= "_{$(c1.site),$(c1.internal_index)}"
        cname1 *= ifelse(!c1.is_annihilation_operator, "^+", "")


        c2 = operator[2]
        cname2 = "C"
        cname2 *= "_{$(c2.site),$(c2.internal_index)}"
        cname2 *= ifelse(!c2.is_annihilation_operator, "^+", "")

        value = h.values[i]
        if value == 1
            valuestring = "+"
        else
            if imag(value) == 0 && real(value) < 0
                valuestring = "-$(-value)"
            else
                valuestring = "+$(value)"
            end
        end

        if value != 0
            cdagcstring *= valuestring * cname1 * cname2 * " "
        end
    end
    println(cdagcstring)
end

using SparseArrays

struct Hamiltonian{T,N} <: AbstractMatrix{T}
    qoperators::QuadraticOPs{T}
    isSC::Bool
    num_internal_degree::Int64
    num_sites::Int64
end

Base.size(h::Hamiltonian{T,N}) where {T,N} = N

function Hamiltonian(num_sites; num_internal_degree=1, isSC=false)
    Hamiltonian(Float64, num_sites; num_internal_degree, isSC)
end

function Hamiltonian(T::DataType, num_sites; num_internal_degree=1, isSC=false)
    qoperators = QuadraticOPs(T)
    N = num_sites * num_internal_degree
    return Hamiltonian{T,N}(qoperators, isSC, num_internal_degree, num_sites)
end

function Base.:+(h::Hamiltonian{T1,N}, term::QuadraticOPs{T2}) where {T1,T2,N}
    qoperators = h.qoperators + term
    Hamiltonian{T1,N}(qoperators, h.isSC, h.num_internal_degree, h.num_sites)
end

using LinearAlgebra

function Base.:*(A::Hamiltonian{T,N}, x::AbstractVector) where {T,N}
    y = zero(x)
    mul!(y, A, x)
    return y
end

#mul!(Y, A, B) -> Y
function LinearAlgebra.mul!(y::AbstractVector, A::Hamiltonian{T,N}, x::AbstractVector) where {T,N}
    y .= 0
    #N = A.num_internal_degree * A.num_sites
    h = A.qoperators
    @inbounds for (i, operator) in enumerate(h.operators)
        c1 = operator[1]
        c2 = operator[2]
        ii = (c1.site - 1) * A.num_internal_degree + c1.internal_index
        jj = (c2.site - 1) * A.num_internal_degree + c2.internal_index
        #ii += ifelse(c1.is_annihilation_operator, N, 0)
        #jj += ifelse(c2.is_annihilation_operator, 0, N)
        y[ii] += h.values[i] * x[jj]
    end
end

#  mul!(C, A, B, α, β) -> C
#A B α + C β
function LinearAlgebra.mul!(y::AbstractVector, A::Hamiltonian{T,N}, x::AbstractVector, α, β) where {T,N}
    y .*= β
    #N = A.num_internal_degree * A.num_sites
    h = A.qoperators
    @inbounds for (i, operator) in enumerate(h.operators)
        c1 = operator[1]
        c2 = operator[2]
        ii = (c1.site - 1) * A.num_internal_degree + c1.internal_index
        jj = (c2.site - 1) * A.num_internal_degree + c2.internal_index
        #ii += ifelse(c1.is_annihilation_operator, N, 0)
        #jj += ifelse(c2.is_annihilation_operator, 0, N)
        y[ii] += α * h.values[i] * x[jj]
    end
end

function Base.display(h::Hamiltonian)
    println("---------------------------------")
    println("Hamiltonian: ")
    println("Num. of sites: $(h.num_sites)")
    println("Num. of internal degree of freedom: $(h.num_internal_degree)")
    if h.isSC
        println("Superconducting state")
    end
    print("H = ")
    display(h.qoperators)
    println("---------------------------------")
end

function check(ham::Hamiltonian{T,N}) where {T,N}
    for (i, operator) in enumerate(h.operators)
        c1 = operator[1]
        c2 = operator[2]
        ii = (c1.site - 1) * ham.num_internal_degree + c1.internal_index
        jj = (c2.site - 1) * ham.num_internal_degree + c2.internal_index
        @assert !c1.is_annihilation_operator && c2.is_annihilation_operator "$i -th term is not C^+ C form: $c1 $c2"
    end
end

function construct_matrix(ham::Hamiltonian{T,N}) where {T,N}
    h = ham.qoperators
    if ham.isSC
        ham_matrix = spzeros(T, 2N, 2N)
    else
        ham_matrix = spzeros(T, N, N)
    end

    for (i, operator) in enumerate(h.operators)
        c1 = operator[1]
        c2 = operator[2]
        ii = (c1.site - 1) * ham.num_internal_degree + c1.internal_index
        jj = (c2.site - 1) * ham.num_internal_degree + c2.internal_index

        #if ham.isSC
        ii += ifelse(c1.is_annihilation_operator, N, 0)
        jj += ifelse(c2.is_annihilation_operator, 0, N)
        #    ham_matrix[ii, jj] += h.values[i]
        #else
        #@assert !c1.is_annihilation_operator && c2.is_annihilation_operator "This is not C^+ C form $c1 $c2"
        ham_matrix[ii, jj] += h.values[i]
        #end
    end
    return ham_matrix
end


struct SCgap{T,N} <: AbstractMatrix{T}
    qoperators::QuadraticOPs{T}
    num_internal_degree::Int64
    num_sites::Int64
end

function SCgap(num_sites; num_internal_degree=1)
    SCgap(Float64, num_sites; num_internal_degree)
end

function SCgap(T::DataType, num_sites; num_internal_degree=1)
    qoperators = QuadraticOPs(T)
    N = num_internal_degree * num_sites
    return SCgap{T,N}(qoperators, num_internal_degree, num_sites)
end

function Base.:+(h::SCgap{T1,N}, term::QuadraticOPs{T2}) where {T1,T2,N}
    qoperators = h.qoperators + term
    SCgap{T1,N}(qoperators, h.num_internal_degree, h.num_sites)
end

function Base.display(h::SCgap)
    println("---------------------------------")
    println("Supercondcting gap: ")
    println("Num. of sites: $(h.num_sites)")
    println("Num. of internal degree of freedom: $(h.num_internal_degree)")
    print("Delta = ")
    display(h.qoperators)
    println("---------------------------------")
end


function construct_matrix(ham::SCgap{T,N}) where {T,N}
    #N = ham.num_internal_degree * ham.num_sites
    h = ham.qoperators

    ham_matrix = spzeros(T, N, N)

    for (i, operator) in enumerate(h.operators)
        c1 = operator[1]
        c2 = operator[2]
        ii = (c1.site - 1) * ham.num_internal_degree + c1.internal_index
        jj = (c2.site - 1) * ham.num_internal_degree + c2.internal_index
        @assert !c1.is_annihilation_operator && !c2.is_annihilation_operator "This is not C^+ C^+ form $c1 $c2"
        ham_matrix[ii, jj] += h.values[i]
    end
    return ham_matrix
end



end
