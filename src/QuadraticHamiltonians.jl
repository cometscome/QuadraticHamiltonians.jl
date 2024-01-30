module QuadraticHamiltonians
export QuadraticOPs, FermionOP, Hamiltonian,
    construct_matrix

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

function Base.display(c::FermionOP{T}) where {T}
    value = c.value
    valuestring = make_valueheader(value)

    cname2 = valuestring * "C"
    cname2 *= "_{$(c.site),$(c.internal_index)}"
    cname2 *= ifelse(!c.is_annihilation_operator, "^+", "")
    println(cname2)
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

function Base.:*(a::T1, h1::QuadraticOPs{T2}) where {T1,T2}
    T = T1
    if T2 <: Complex
        T = ComplexF64
    end
    operators = Tuple{FermionOP{T},FermionOP{T}}[]
    values = T[]
    for (i, operator) in enumerate(h1.operators)
        c1 = operator[1]
        c2 = operator[2]
        value = a * h1.values[i]
        push!(operators, (FermionOP(c1.site, c1.internal_index, c1.is_annihilation_operator, one(T)),
            FermionOP(c2.site, c2.internal_index, c2.is_annihilation_operator, one(T))))
        push!(values, value)
    end
    return QuadraticOPs(operators, values)
end

function Base.:*(h1::QuadraticOPs{T1}, a::Number) where {T1}
    return Base.:*(a::Number, h1::QuadraticOPs{T1})
end

function Base.:+(h1::QuadraticOPs{T1}, h2::QuadraticOPs{T2}) where {T1,T2}
    return quadratic_plus(h1::QuadraticOPs{T1}, h2::QuadraticOPs{T2}, +1)
end

function Base.:-(h1::QuadraticOPs{T1}, h2::QuadraticOPs{T2}) where {T1,T2}
    return quadratic_plus(h1::QuadraticOPs{T1}, h2::QuadraticOPs{T2}, -1)
end

function quadratic_plus(h1::QuadraticOPs{T1}, h2::QuadraticOPs{T2}, sign::Number) where {T1,T2}
    T = T1
    if T2 <: Complex
        T = ComplexF64
    end
    if T <: Int && T2 <: Real
        T = T2
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

        value = sign * h2.values[i]

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


function make_valueheader(value)

    if value == 1
        valuestring = "+"
    else
        if imag(value) == 0
            if real(value) < 0
                valuestring = "-$(-real(value))"
            else
                valuestring = "+$(real(value))"
            end
        elseif real(value) == 0
            if imag(value) < 0
                valuestring = "-$(-imag(value))im"
            else
                valuestring = "+$(imag(value))im"
            end
        else
            if imag(value) != 0
                valuestring = "+($(value))"
            else
                valuestring = "+$(real(value))"
            end
        end
    end
    return valuestring
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
        valuestring = make_valueheader(value)



        if value != 0
            cdagcstring *= valuestring * cname1 * cname2 * " "
        end
    end
    println(cdagcstring)
end

using SparseArrays

struct Hamiltonian{T,N,isSC,num_internal_degree,num_sites} <: AbstractMatrix{T}
    matrix::SparseMatrixCSC{T,Int64}
end

Base.size(h::Hamiltonian) = size(h.matrix)
Base.getindex(A::Hamiltonian, i::Int) = getindex(A.matrix, i)
Base.getindex(A::Hamiltonian, I::Vararg{Int,N}) where {N} = getindex(A.matrix, I)


function Hamiltonian(num_sites; num_internal_degree=1, isSC=false)
    Hamiltonian(Float64, num_sites; num_internal_degree, isSC)
end

function Hamiltonian(T::DataType, num_sites; num_internal_degree=1, isSC=false)
    N = num_sites * num_internal_degree
    if isSC
        matrix = spzeros(T, 2N, 2N)
    else
        matrix = spzeros(T, N, N)
    end
    return Hamiltonian{T,N,isSC,num_internal_degree,num_sites}(matrix)
end

function hamiltonian_plus(h::Hamiltonian{T1,N,isSC,num_internal_degree,num_sites}, term::QuadraticOPs{T2}, sign::Number) where {T1,T2,N,isSC,num_internal_degree,num_sites}
    for (i, operator) in enumerate(term.operators)
        c1 = operator[1]
        c2 = operator[2]
        @assert c1.internal_index <= num_internal_degree "internal degree of freedom is $num_internal_degree but the internal index of the $i -th first operator is $(c1.internal_index)"
        @assert c2.internal_index <= num_internal_degree "internal degree of freedom is $num_internal_degree but the internal index of the $i -th second operator is $(c2.internal_index)"

        if isSC
            ii = (c1.site - 1) * num_internal_degree + c1.internal_index + c1.is_annihilation_operator * N
            jj = (c2.site - 1) * num_internal_degree + c2.internal_index + (!c2.is_annihilation_operator) * N
        else
            @assert !c1.is_annihilation_operator && c2.is_annihilation_operator "$i -th term is not C^+ C form $c1 $c2"
            ii = (c1.site - 1) * num_internal_degree + c1.internal_index
            jj = (c2.site - 1) * num_internal_degree + c2.internal_index
        end
        h.matrix[ii, jj] += sign * term.values[i]
    end
    Hamiltonian{T1,N,isSC,num_internal_degree,num_sites}(h.matrix)
end

function Base.:+(h::Hamiltonian{T1,N,isSC,num_internal_degree,num_sites}, term::QuadraticOPs{T2}) where {T1,T2,N,isSC,num_internal_degree,num_sites}
    return hamiltonian_plus(h, term, +1)
end

function Base.:-(h::Hamiltonian{T1,N,isSC,num_internal_degree,num_sites}, term::QuadraticOPs{T2}) where {T1,T2,N,isSC,num_internal_degree,num_sites}
    return hamiltonian_plus(h, term, -1)
end

function get_coefficient(h::Hamiltonian{T1,N,isSC,num_internal_degree,num_sites},
    c1::FermionOP, c2::FermionOP) where {T1,N,isSC,num_internal_degree,num_sites}

    ii = (c1.site - 1) * num_internal_degree + c1.internal_index
    jj = (c2.site - 1) * num_internal_degree + c2.internal_index
    if isSC
        ii += ifelse(c1.is_annihilation_operator, N, 0)
        jj += ifelse(c2.is_annihilation_operator, 0, N)
    else
        @assert !c1.is_annihilation_operator && c2.is_annihilation_operator "This is not C^+ C form $c1 $c2"
    end
    return h.matrix[ii, jj]
end

Base.getindex(A::Hamiltonian, c1::FermionOP, c2::FermionOP) = get_coefficient(A, c1, c2)



using LinearAlgebra

function Base.:*(A::Hamiltonian, x::AbstractVector)
    y = zero(x)
    mul!(y, A, x)
    return y
end

#mul!(Y, A, B) -> Y
function LinearAlgebra.mul!(y::AbstractVector, A::Hamiltonian{T,N,true,num_internal_degree,num_sites}, x::AbstractVector) where {T,N,num_internal_degree,num_sites}
    mul!(y, A.matrix, x)
end

function LinearAlgebra.mul!(y::AbstractVector, A::Hamiltonian{T,N,false,num_internal_degree,num_sites}, x::AbstractVector) where {T,N,num_internal_degree,num_sites}
    mul!(y, A.matrix, x)
end

#  mul!(C, A, B, α, β) -> C
#A B α + C β
function LinearAlgebra.mul!(y::AbstractVector, A::Hamiltonian{T,N,true,num_internal_degree,num_sites}, x::AbstractVector, α::Number, β::Number) where {T,N,num_internal_degree,num_sites}
    mul!(y, A.matrix, x, α, β)

end


function Base.display(h::Hamiltonian{T,N,isSC,num_internal_degree,num_sites}) where {T,N,num_internal_degree,num_sites,isSC}
    println("---------------------------------")
    println("Hamiltonian: ")
    println("Num. of sites: $(num_sites)")
    println("Num. of internal degree of freedom: $(num_internal_degree)")
    if isSC
        println("Superconducting state")
    end
    print("H = ")
    display_quadratic(h)
    #display(h.qoperators)
    println("---------------------------------")
end

function display_quadratic(ham::Hamiltonian{T,N,true,num_internal_degree,num_sites}) where {T,N,num_internal_degree,num_sites}
    cdagcstring = ""
    for j = 1:2N
        if j <= N
            jinternal = (j - 1) % num_internal_degree + 1
            jsite = (j - jinternal) ÷ num_internal_degree + 1
            j_isdag = false
        else
            jinternal = (j - N - 1) % num_internal_degree + 1
            jsite = (j - N - jinternal) ÷ num_internal_degree + 1
            j_isdag = true
        end

        for ii = ham.matrix.colptr[j]:ham.matrix.colptr[j+1]-1

            if ii <= length(ham.matrix.rowval)
                i = ham.matrix.rowval[ii]
                if i <= N
                    iinternal = (i - 1) % num_internal_degree + 1
                    isite = (i - iinternal) ÷ num_internal_degree + 1
                    i_isdag = true
                else
                    iinternal = (i - N - 1) % num_internal_degree + 1
                    isite = (i - N - iinternal) ÷ num_internal_degree + 1
                    i_isdag = false
                end


                cname1 = "C"
                cname1 *= "_{$(isite),$(iinternal)}"
                cname1 *= ifelse(i_isdag, "^+", "")


                cname2 = "C"
                cname2 *= "_{$(jsite),$(jinternal)}"
                cname2 *= ifelse(j_isdag, "^+", "")

                value = ham.matrix.nzval[ii]
                valuestring = make_valueheader(value)


                if value != 0
                    cdagcstring *= valuestring * cname1 * cname2 * " "
                end
            end
        end


    end
    println(cdagcstring)
end


function display_quadratic(ham::Hamiltonian{T,N,false,num_internal_degree,num_sites}) where {T,N,num_internal_degree,num_sites}
    cdagcstring = ""
    for j = 1:N
        jinternal = (j - 1) % num_internal_degree + 1
        jsite = (j - jinternal) ÷ num_internal_degree + 1
        j_isdag = false


        for ii = ham.matrix.colptr[j]:ham.matrix.colptr[j+1]-1
            if ii <= length(ham.matrix.rowval)
                i = ham.matrix.rowval[ii]

                iinternal = (i - 1) % num_internal_degree + 1
                isite = (i - iinternal) ÷ num_internal_degree + 1
                i_isdag = true


                cname1 = "C"
                cname1 *= "_{$(isite),$(iinternal)}"
                cname1 *= ifelse(i_isdag, "^+", "")


                cname2 = "C"
                cname2 *= "_{$(jsite),$(jinternal)}"
                cname2 *= ifelse(j_isdag, "^+", "")

                value = ham.matrix.nzval[ii]
                valuestring = make_valueheader(value)

                if value != 0
                    cdagcstring *= valuestring * cname1 * cname2 * " "
                end
            end
        end


    end
    println(cdagcstring)
end


function construct_matrix(ham)
    return ham.matrix
end


include("RSCGSolver.jl")
export RSCGSolver, solve

include("ChebyshevSolver.jl")
export ChebyshevSolver

include("Meanfields.jl")
export Meanfields_solver

end
