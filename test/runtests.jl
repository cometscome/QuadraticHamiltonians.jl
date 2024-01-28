using QuadraticHamiltonians
using Test
using LinearAlgebra
using BenchmarkTools
using KrylovKit
using SparseArrays
using InteractiveUtils

function test4()
    μ = -1.5
    N = 16
    Δ = 0.5
    ham = Hamiltonian(N, num_internal_degree=2, isSC=true)
    hops = [+1, -1]
    Δ = 0.5
    for i = 1:N
        for ispin = 1:2
            ci = FermionOP(i, ispin)
            for d in hops
                j = i + d
                j += ifelse(j > N, -N, 0)
                j += ifelse(j < 1, N, 0)
                jspin = ispin
                cj = FermionOP(j, jspin)
                ham += -1 * (ci' * cj - ci * cj')
            end
            ham += -μ * (ci' * ci - ci * ci')
        end
        ciup = FermionOP(i, 1)
        cidown = FermionOP(i, 2)
        ham += Δ * ciup' * cidown' + Δ * cidown * ciup
        ham += -Δ * cidown' * ciup' - Δ * ciup * cidown
    end
    T = 0.1
    m = Meanfields_solver(ham, T)
    c1up = FermionOP(1, 1)
    c1down = FermionOP(1, 2)
    Gij0 = calc_meanfields(m, c1up, c1down) #<c1up c1down>
    println(Gij0)
    num = 100
    zs = zeros(ComplexF64, num)
    ene = range(-1, 1, length=num)
    eta = 0.01
    for i = 1:num
        zs[i] = ene[i] + im * eta
    end
    Gij = calc_greenfunction(ham, zs, c1up', c1up)
    fp = open("test.txt", "w")
    for i = 1:num
        #println(Gij[i])
        println(fp, ene[i], "\t", -imag(Gij[i]) / pi)
    end
    close(fp)
end

function test3()
    μ = -1.5
    N = 16
    Δ = 0.5
    ham = Hamiltonian(N, isSC=true)
    hops = [+1, -1]
    Δ = 0.5
    for i = 1:N
        ci = FermionOP(i)
        for d in hops
            j = i + d
            j += ifelse(j > N, -N, 0)
            j += ifelse(j < 1, N, 0)
            cj = FermionOP(j)
            ham += -1 * (ci' * cj - ci * cj')
        end
        ham += -μ * (ci' * ci - ci * ci')
        ham += Δ * ci' * ci' + Δ * ci * ci
    end
    T = 0.1
    m = Meanfields_solver(ham, T)
    c1 = FermionOP(1)
    Gij0 = calc_meanfields(m, c1, c1) #<c1 c1>
    println(Gij0)
end

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
    #test2()
    #test3()
    test4()
end
