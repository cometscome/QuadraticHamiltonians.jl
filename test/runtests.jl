using QuadraticHamiltonians
using Test

function test2()
    N = 16
    ham = Hamiltonian(N)
    t = -1
    for i = 1:N
        ci = FermionOP(i)
        j = i + 1
        j += ifelse(j > N, -N, 0)
        cj = FermionOP(j)
        ham += t * ci' * cj

        j = i - 1
        j += ifelse(j < 1, N, 0)

        cj = FermionOP(j)
        ham += t * ci' * cj
    end
    display(ham)
    ham_matrix = construct_matrix(ham)
    display(ham_matrix)
end

function test()
    c1 = FermionOP(1)
    c2 = FermionOP(2)

    N = 16
    ham = Hamiltonian(N)
    println(ham)
    ham += 10 * c1' * c1 + 3 * c2' * c2 + 3 * c1' * c2
    display(ham)
    ham_matrix = construct_matrix(ham)
    display(ham_matrix)
    return

    c1 = FermionOP(1)
    c1dag = c1'

    c2 = FermionOP(2)
    c1dag = c2'
    h = 2 * c1' * c1 + 3 * c2' * c2 + 3 * c1' * c2
    hd = h'
    display(h)
    display(hd)

    println(QuadraticHamiltonians.find_operator(h, (c2', c2)))

    h += 4 * c1' * c2
    display(h)
end

@testset "QuadraticHamiltonians.jl" begin
    # Write your tests here.
    #test()
    test2()
end