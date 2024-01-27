using QuadraticHamiltonians
using Test

function test()
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
    test()
end
