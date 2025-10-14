
using SparseArrays

struct Hamiltonian{T,N,isSC,num_internal_degree,num_sites} <: AbstractMatrix{T}
    matrix::SparseMatrixCSC{T,Int64}
end

Base.size(h::Hamiltonian) = size(h.matrix)
Base.getindex(A::Hamiltonian, i::Int) = getindex(A.matrix, i)
Base.getindex(A::Hamiltonian, I::Vararg{Int,N}) where {N} = getindex(A.matrix, I)
Base.getindex(A::Hamiltonian, I...) = getindex(A.matrix, I...)


function Base.setindex!(A::Hamiltonian, v, i::Int)
    A.matrix[i] = v
end
function Base.setindex!(A::Hamiltonian, v, I::Vararg{Int,N}) where {N}
    #println(I)
    A.matrix[I...] = v
end
function Base.setindex!(A::Hamiltonian, X, I...)
    A.matrix[I...] = X
end

function get_dim(A::Hamiltonian)
    n, _ = size(A.matrix)
    return n
end

function Hamiltonian(num_sites; num_internal_degree=1, isSC=false)
    Hamiltonian(Float64, num_sites; num_internal_degree, isSC)
end

function Hamiltonian(matrix::AbstractMatrix, num_sites; num_internal_degree=1)
    matrixSP = sparse(matrix)
    dim, _ = size(matrix)
    if dim == 2 * num_sites * num_internal_degree
        isSC = true
    elseif dim == num_sites * num_internal_degree
        isSC = false
    else
        error("size of matrix is wrong")
    end
    T = eltype(matrix)
    N = num_sites * num_internal_degree
    return Hamiltonian{T,N,isSC,num_internal_degree,num_sites}(matrixSP)
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

function Base.setindex!(h::Hamiltonian{T1,N,isSC,num_internal_degree,num_sites}, v,
    c1::FermionOP, c2::FermionOP) where {T1,N,isSC,num_internal_degree,num_sites}

    ii = (c1.site - 1) * num_internal_degree + c1.internal_index
    jj = (c2.site - 1) * num_internal_degree + c2.internal_index
    if isSC
        ii += ifelse(c1.is_annihilation_operator, N, 0)
        jj += ifelse(c2.is_annihilation_operator, 0, N)
    else
        @assert !c1.is_annihilation_operator && c2.is_annihilation_operator "This is not C^+ C form $c1 $c2"
    end
    h.matrix[ii, jj] = v
end


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

indexer_x_2D(i, Nx) = (i - 1) % Nx + 1
function indexer_y_2D(i, Nx)
    return div(i - 1, Nx) + 1
end

function make_X(ham::Hamiltonian{T,N,isSC,num_internal_degree,num_sites}, Nx) where {T,N,isSC,num_internal_degree,num_sites}
    X = make_X(ham, Nx, indexer_x_2D)
    return X
end

function make_X(ham::Hamiltonian{T,N,isSC,num_internal_degree,num_sites}, Nx, indexer_x) where {T,N,isSC,num_internal_degree,num_sites}
    X = zero(ham.matrix)
    for i = 1:num_sites
        ix = indexer_x(i, Nx)
        for i_i = 1:num_internal_degree
            ii = (i - 1) * num_internal_degree + i_i
            jj = ii
            v = ix
            X[ii, jj] = v
            if isSC
                X[ii+N, jj+N] = v
            end
        end
    end
    return X
end

function make_Y(ham::Hamiltonian{T,N,isSC,num_internal_degree,num_sites}, Nx) where {T,N,isSC,num_internal_degree,num_sites}
    Y = make_Y(ham, Nx, indexer_y_2D)
    return Y
end

function make_Y(ham::Hamiltonian{T,N,isSC,num_internal_degree,num_sites}, Nx, indexer_y) where {T,N,isSC,num_internal_degree,num_sites}
    Y = zero(ham.matrix)
    for i = 1:num_sites
        iy = indexer_y(i, Nx)
        for i_i = 1:num_internal_degree
            ii = (i - 1) * num_internal_degree + i_i
            jj = ii
            v = iy
            Y[ii, jj] = v
            if isSC
                Y[ii+N, jj+N] = v
            end
        end
    end
    return Y
end
