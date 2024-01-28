using QuadraticHamiltonians
using Test
using LinearAlgebra
using BenchmarkTools
using KrylovKit
using SparseArrays
using InteractiveUtils

function test2()
    T = 0.1
    μ = -1.5
    rscg = RSCGSolver(T)



    Nx = 16 * 4
    Ny = 16 * 4
    N = Nx * Ny
    Δ = 0.5
    ham = Hamiltonian(N, isSC=true)
    ham_test = spzeros(2N, 2N)
    t = -1
    for ix = 1:Nx
        for iy = 1:Ny
            i = (iy - 1) * Nx + ix
            ci = FermionOP(i)

            jx = ix + 1
            jx += ifelse(jx > Nx, -Nx, 0)
            jy = iy
            j = (jy - 1) * Nx + jx
            cj = FermionOP(j)
            ham += t * (ci' * cj - ci * cj')
            ham_test[i, j] = t
            ham_test[i+N, j+N] = -t

            jx = ix - 1
            jx += ifelse(jx < 1, Nx, 0)
            jy = iy
            j = (jy - 1) * Nx + jx
            cj = FermionOP(j)
            ham += t * (ci' * cj - ci * cj')
            ham_test[i, j] = t
            ham_test[i+N, j+N] = -t

            jy = iy + 1
            jy += ifelse(jy > Ny, -Ny, 0)
            jx = ix
            j = (jy - 1) * Nx + jx
            cj = FermionOP(j)
            ham += t * (ci' * cj - ci * cj')
            ham_test[i, j] = t
            ham_test[i+N, j+N] = -t

            jy = iy - 1
            jy += ifelse(jy < 1, Ny, 0)
            jx = ix
            j = (jy - 1) * Nx + jx
            cj = FermionOP(j)
            ham += t * (ci' * cj - ci * cj')
            ham_test[i, j] = t
            ham_test[i+N, j+N] = -t

            j = i
            cj = FermionOP(j)
            ham += -μ * (ci' * cj - ci * cj')
            ham_test[i, j] = -μ
            ham_test[i+N, j+N] = μ

            ham += Δ * ci' * cj' + Δ * ci * cj
            ham_test[i, j+N] = Δ
            ham_test[i+N, j] = Δ
        end
    end

    m = Meanfields_solver(ham, T, wmax=12.0)
    i = 1
    j = 1 + N
    ham_matrix = construct_matrix(ham)
    println(ham_matrix - ham_matrix')
    println(ham_matrix - ham_test)
    display(ham_matrix)
    display(ham_test)

    #=
    x = rand(2N)
    y = rand(2N)
    @time mul!(y, ham, x)
    @time mul!(y, ham, x)
    @code_warntype mul!(y, ham, x)
    return
    =#


    println("sparse matrix")
    @btime solve($rscg, $ham_matrix, $i, $j)
    Gij0 = solve(rscg, ham_matrix, i, j)
    println(Gij0)

    #@code_warntype solve(rscg, ham_matrix, i, j)

    println("hamiltonian,operator")
    c1 = FermionOP(1)
    c2 = FermionOP(2)
    @btime calc_meanfields($m, $(c1), $c1)
    Gij0 = calc_meanfields(m, c1, c1)
    println(Gij0)

    #@code_warntype calc_meanfields(m, c1, c1)
    # return

    println("hamiltonian,ij")
    @btime calc_meanfields($m, $i, $j)
    Gij0 = calc_meanfields(m, i, j)
    println(Gij0)
    #@btime calc_meanfields($m, $(c1'), $c1)
    #println(Gij0)
    return

    x = rand(N)
    vals, vecs, info = eigsolve(ham, x, 1, :SR, ishermitian=true)
    println(vals)
    return

    display(ham)
    ham_matrix = construct_matrix(ham)
    display(ham_matrix)

    x = rand(N)
    y = zero(x)

    @btime mul!($y, $ham, $x)
    println(sum(x))

    println("matrix")
    @btime mul!($y, $ham_matrix, $x)
    println(sum(x))

end

function test()
    c1 = FermionOP(1)
    c2 = FermionOP(2)

    N = 16
    ham = Hamiltonian(N)
    println(ham)
    ham += 10 * c1' * c1 + 3 * c2' * c2 + 3 * c1' * c2
    display(ham)
    return
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
