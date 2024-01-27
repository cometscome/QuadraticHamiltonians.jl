module QuadraticHamiltonians
export QuadraticOPs, FermionOP
# Write your package code here.


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
    return FermionOP(site, internal_index, false, 1)
end
function FermionOP(site)
    return FermionOP(site, 1, false, 1)
end

Base.adjoint(c::FermionOP) = FermionOP(c.site, c.internal_index, !c.is_annihilation_operator, c.value)

struct QuadraticOPs{T}
    operators::Vector{Tuple{FermionOP{T},FermionOP{T}}}
    values::Vector{T}
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
        push!(operators, (FermionOP(c1.site, c1.internal_index, c1.is_annihilation_operator, one(value)),
            FermionOP(c2.site, c2.internal_index, c2.is_annihilation_operator, one(value))))
        push!(values, value)
    end

    for (i, operator) in enumerate(h2.operators)
        c1 = operator[1]
        c2 = operator[2]
        value = h2.values[i]
        op = (FermionOP(c1.site, c1.internal_index, c1.is_annihilation_operator, one(value)),
            FermionOP(c2.site, c2.internal_index, c2.is_annihilation_operator, one(value)))
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
        cname1 *= ifelse(c1.is_annihilation_operator, "^+", "")


        c2 = operator[2]
        cname2 = "C"
        cname2 *= "_{$(c2.site),$(c2.internal_index)}"
        cname2 *= ifelse(c2.is_annihilation_operator, "^+", "")

        value = h.values[i]
        if value == 1
            valuestring = "+"
        else
            valuestring = "+$(value)"
        end

        if value != 0
            cdagcstring *= valuestring * cname1 * cname2
        end
    end
    println(cdagcstring)
end


end
