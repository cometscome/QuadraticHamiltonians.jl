module LK
using LinearAlgebra
using SparseArrays

mutable struct LKvector{T<:Number} <: AbstractVector{T}
    values::Array{T,1}
    N::Int64
    indices::Array{Int64,1}
    nonzeros::Int64
    hasvalue::Vector{Bool}
    eps::Float64
end
export LKvector, clear!, filter!

Base.size(A::LKvector) = (A.N,)
Base.length(A::LKvector) = A.N

function LKvector(B::Type, N, eps)
    values = zeros(B, N)
    indices = zeros(Int64, N)
    hasvalue = zeros(Bool, N)
    nonzeros = 0
    return LKvector{B}(values, N, indices, nonzeros, hasvalue, eps)
end

function LKvector(B::Type, N)
    values = zeros(B, N)
    indices = zeros(Int64, N)
    hasvalue = zeros(Bool, N)
    nonzeros = 0
    eps = 1e-16
    return LKvector{B}(values, N, indices, nonzeros, hasvalue, eps)
end

function LKvector(N)
    return LKvector(ComplexF64, N, 1e-16)
end

function LKvector(N, eps)
    return LKvector(ComplexF64, N, eps)
end

function LinearAlgebra.norm(A::LKvector)
    v = 0
    @inbounds for ip = 1:A.nonzeros
        i = A.indices[ip]
        v += abs(A.values[i])^2
        #conj(A.values[i])*A.values[i]
    end
    return sqrt(v)
end

function LinearAlgebra.dot(a::LKvector, b::LKvector)
    nonzeros = a.nonzeros
    aorb = 1
    v = 0
    if nonzeros < b.nonzeros
        nonzeros = b.nonzeros
        aorb = 2
    end

    for ip = 1:nonzeros
        if aorb == 1
            i = a.indices[ip]
        else
            i = b.indices[ip]
        end
        ai = a.values[i]
        bi = b.values[i]
        v += conj(ai) * bi
    end
    return v
end




function Base.print(A::LKvector)
    for ip = 1:A.nonzeros
        i = A.indices[ip]
        println(i, "\t", A.values[i])
    end
end

Base.zeros(::Type{LKvector}, N::Union{Integer,AbstractUnitRange}) = LKvector(N)
Base.zero(A::LKvector) = LKvector(eltype(A.values), A.N, A.eps)

function Base.copy(A::LKvector)
    B = LKvector(eltype(A.values), A.N, A.eps)
    B.nonzeros = A.nonzeros
    for ip = 1:A.nonzeros
        i = A.indices[ip]
        B.indices[ip] = i
        B.values[i] = A.values[i]
    end
    return B
end

function Base.copy!(B::LKvector, A::LKvector)
    #B = LKvector(eltype(A.values),A.N,A.eps)
    B.nonzeros = A.nonzeros
    for ip = 1:A.nonzeros
        i = A.indices[ip]
        B.indices[ip] = i
        B.values[i] = A.values[i]
    end
    #return B
end


Base.getindex(A::LKvector, i::Int) = A.values[i]

function Base.setindex!(A::LKvector, v, i::Int)
    #println("val",A.values[i])
    #println("v",v)
    if A.hasvalue[i] == false
        #if abs(A.values[i]) <= 0
        if abs(v) < A.eps
            return
        else
            A.nonzeros += 1
            A.indices[A.nonzeros] = i
            A.values[i] = v
            A.hasvalue[i] = true
            #println("nonzeros, ",A.nonzeros)
            #println("index, ",A.indices[A.nonzeros] )
        end
    else
        if abs(v) < A.eps
            A.values[i] = 0
            A.hasvalue[i] = false
            temp = A.indices[1:A.nonzeros]
            filter!(x -> x != i, temp)
            A.indices[A.nonzeros] = 0
            A.nonzeros -= 1
            A.indices[1:A.nonzeros] = temp[1:A.nonzeros]
            return
        else
            A.values[i] = v
        end
    end

end

function Base.:*(a::Number, x::LKvector)
    y = copy(x)
    for ip = 1:x.nonzeros
        i = x.indices[ip]
        y.values[i] = a * x.values[i]
    end
    return y
end

function LinearAlgebra.mul!(a::Number, x::LKvector)
    for ip = 1:x.nonzeros
        i = x.indices[ip]
        x.values[i] = a * x.values[i]
    end
end

function Base.:*(x::LKvector, a::Number)
    y = copy(x)
    for ip = 1:x.nonzeros
        i = x.indices[ip]
        y.values[i] = a * x.values[i]
    end
    return y
end

function Base.:*(A::SparseArrays.SparseMatrixCSC, x::LKvector)
    y = LKvector(x)
    for ip = 1:x.nonzeros #y_i = sum_l A_{il} x_l
        l = x.indices[ip]
        xl = x.values[l]

        for ii = A.colptr[l]:(A.colptr[l+1]-1)
            i = A.rowval[ii]
            Ail = A.nzval[ii]
            add_i!(y, i, Ail * xl)
        end
    end
    return y
end

function LinearAlgebra.mul!(C::LKvector, A::SparseArrays.SparseMatrixCSC, B::LKvector, α::Number, β::Number) #A B α + C β
    #println("C1 ", norm(C))
    mul!(β, C)
    #println("C2 ", norm(C))
    @inbounds for ip = 1:B.nonzeros
        l = B.indices[ip]
        if B.hasvalue[l]
            xl = B.values[l]
            #if abs(xl) != 0
            for ii = A.colptr[l]:(A.colptr[l+1]-1)
                i = A.rowval[ii]
                Ail = A.nzval[ii]
                add_i!(C, i, α * Ail * xl)
            end
            #end
        end
        #println("C3 ", norm(C))
    end
    #r = rand()
    #if r < 0.5
    filter!(C)
    #end

end

function LinearAlgebra.mul!(C::LKvector, A::SparseArrays.SparseMatrixCSC, B::LKvector) #C = A B α
    #println("C1 ", norm(C))
    #mul!(β,C) 
    clear!(C)
    #println("C2 ", norm(C))
    for ip = 1:B.nonzeros
        l = B.indices[ip]
        if B.hasvalue[l]
            xl = B.values[l]
            #if abs(xl) != 0
            for ii = A.colptr[l]:(A.colptr[l+1]-1)
                i = A.rowval[ii]
                Ail = A.nzval[ii]
                add_i!(C, i, Ail * xl)
            end
            #end
        end
        #println("C3 ", norm(C))
    end
    #r = rand()
    #if r < 0.5
    #filter!(C)
    #end
    #println("done") 
end

function clear!(A::LKvector)
    #println("A.nonzeros ",A.nonzeros)
    for ip = 1:A.nonzeros
        i = A.indices[ip]
        A.values[i] = 0
        A.indices[ip] = 0
    end
    A.nonzeros = 0
end

function Base.filter!(A::LKvector)
    #println("A.nonzeros ",A.nonzeros)
    count = 0
    #indices = zeros(Int64, A.nonzeros)
    #println("before ", A.nonzeros)

    for ip = 1:A.nonzeros
        i = A.indices[ip]
        #A.indices[ip] = 0
        v = A.values[i]
        if abs(v) <= A.eps
            A.values[i] = 1e-16
            A.hasvalue[i] = false
            #println("filtered")
        else
            count += 1
            #indices[count] = i
            A.indices[count] = i
        end
    end
    #println(count)
    A.indices[count+1:A.nonzeros] .= 0

    #for i = 1:count
    #    A.indices[i] = indices[i]
    #end

    A.nonzeros = count
    #println("after ", A.nonzeros)
end


@inline function add_i!(x::LKvector, i, v)
    #println(i,"\t",x.nonzeros,"\t",abs(x.values[i]))
    if abs(v) > x.eps
        if x.hasvalue[i] == false
            #if abs(x.values[i]) <= 1e-16
            x.nonzeros += 1
            #println(x.nonzeros, "\t", x.values[i], "\t", abs(v))
            x.indices[x.nonzeros] = i
            x.values[i] = v #+ 1e-16
            x.hasvalue[i] = true
        else
            x.values[i] += v
            #if abs(x.values[i]) <= 0
            #    println("near zero, ", x.values[i])
            #end
            #=
            if abs(x.values[i]) < x.eps
                x.values[i] = 0
                x.hasvalue[i] = false
                ip = findfirst(x -> x == i, x.indices)
                #temp = x.indices[1:x.nonzeros]
                #filter!(x -> x != i, temp)
                for k = ip:x.nonzeros-1
                    x.indices[k] = x.indices[k+1]
                end
                x.indices[x.nonzeros] = 0
                x.nonzeros -= 1

                #x.indices[1:x.nonzeros] = temp[1:x.nonzeros]
            end
            =#

        end
    end
end


function LKvector(x::LKvector)
    return LKvector(eltype(x.values), x.N, x.eps)
end

function LKvector(x::AbstractVector)
    eps = 1e-16
    return LKvector(x, eps)
end

function LKvector(x::AbstractVector, eps)
    N = length(x)
    xout = LKvector(eltype(x), N, eps)
    for i = 1:N
        if abs(x[i]) >= eps
            xout.nonzeros += 1
            xout.values[i] = x[i]
            xout.indices[xout.nonzeros] = i
        end
    end
    return xout
end





end