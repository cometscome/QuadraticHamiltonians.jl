using QuadraticHamiltonians
using Test
using LinearAlgebra
using BenchmarkTools
using KrylovKit
using SparseArrays
using InteractiveUtils

const σ0 = [1 0
    0 1]
const σx = [0 1
    1 0]
const σy = [0 -im
    im 0]
const σz = [1 0
    0 -1]

function make_x_plushop(Nx, Ny, BC)
    N = Nx * Ny
    Txhop = spzeros(Int64, N, N)
    for ix = 1:Nx
        for iy = 1:Ny
            i = (iy - 1) * Nx + ix
            jx = ix + 1
            jy = iy
            if BC == "PBC"
                jx += ifelse(jx > Nx, -Nx, 0)
            elseif BC == "OBC"
            else
                error("BC = $BC is not supported. Use PBC or OBC")
            end
            if 1 <= jx <= Nx
                j = (jy - 1) * Nx + jx
                Txhop[i, j] = 1
            end
        end
    end
    return Txhop
end

function make_x_minushop(Nx, Ny, BC)
    N = Nx * Ny
    Txhop = spzeros(Int64, N, N)
    for ix = 1:Nx
        for iy = 1:Ny
            i = (iy - 1) * Nx + ix
            jx = ix - 1
            jy = iy
            if BC == "PBC"
                jx += ifelse(jx < 1, Nx, 0)
            elseif BC == "OBC"
            else
                error("BC = $BC is not supported. Use PBC or OBC")
            end
            if 1 <= jx <= Nx
                j = (jy - 1) * Nx + jx
                Txhop[i, j] = 1
            end
        end
    end
    return Txhop
end

function make_y_plushop(Nx, Ny, BC)
    N = Nx * Ny
    Tyhop = spzeros(Int64, N, N)
    for ix = 1:Nx
        for iy = 1:Ny
            i = (iy - 1) * Nx + ix
            jx = ix
            jy = iy + 1
            if BC == "PBC"
                jy += ifelse(jy > Ny, -Ny, 0)
            elseif BC == "OBC"
            else
                error("BC = $BC is not supported. Use PBC or OBC")
            end
            if 1 <= jy <= Ny
                j = (jy - 1) * Nx + jx
                Tyhop[i, j] = 1
            end
        end
    end
    return Tyhop
end

function make_y_minushop(Nx, Ny, BC)
    N = Nx * Ny
    Tyhop = spzeros(Int64, N, N)
    for ix = 1:Nx
        for iy = 1:Ny
            i = (iy - 1) * Nx + ix
            jx = ix
            jy = iy - 1
            if BC == "PBC"
                jy += ifelse(jy < 1, Ny, 0)
            elseif BC == "OBC"
            else
                error("BC = $BC is not supported. Use PBC or OBC")
            end
            if 1 <= jy <= Ny
                j = (jy - 1) * Nx + jx
                Tyhop[i, j] = 1
            end
        end
    end
    return Tyhop
end


function make_TSC_hamiltonian(Nx, Ny, μ, Δs, h, α, isxPBC, isyPBC)
    N = Nx * Ny
    ham = Hamiltonian(ComplexF64, N, num_internal_degree=2, isSC=true)
    hops = [(+1, 0), (-1, 0), (0, +1), (0, -1)]
    for ix = 1:Nx
        for iy = 1:Ny
            i = (iy - 1) * Nx + ix
            for ispin = 1:2
                ci = FermionOP(i, ispin)
                σy = ifelse(ispin == 1, -im, im)
                σx = 1
                σz = ifelse(ispin == 1, 1, -1)

                for (dx, dy) in hops
                    jx = ix + dx
                    if isxPBC
                        jx += ifelse(jx > Nx, -Nx, 0)
                        jx += ifelse(jx < 1, Nx, 0)
                    end
                    jy = iy + dy
                    if isyPBC
                        jy += ifelse(jy > Ny, -Ny, 0)
                        jy += ifelse(jy < 1, Ny, 0)
                    end

                    if 1 <= jx <= Nx && 1 <= jy <= Ny
                        j = (jy - 1) * Nx + jx
                        jspin = ispin
                        cj = FermionOP(j, jspin)

                        ham += -1 * (ci' * cj - cj * ci')

                        jspin = ifelse(ispin == 1, 2, 1)
                        cj = FermionOP(j, jspin)

                        if dy == 0
                            if dx == 1
                                ham += (α / (2im)) * (ci' * cj - cj * ci') * σy
                            elseif dx == -1
                                ham += -1 * (α / (2im)) * (ci' * cj - cj * ci') * σy
                            end
                        elseif dx == 0
                            if dy == 1
                                ham += (α / (2im)) * (ci' * cj - cj * ci') * σx
                            elseif dy == -1
                                ham += -1 * (α / (2im)) * (ci' * cj - cj * ci') * σx
                            end
                        end
                    end
                end
                ham += (-μ - h * σz) * (ci' * ci - ci * ci')
            end
            ciup = FermionOP(i, 1)
            cidown = FermionOP(i, 2)
            ham += Δs[i] * ciup' * cidown' + Δs[i] * cidown * ciup
            ham += -Δs[i] * cidown' * ciup' - Δs[i] * ciup * cidown
        end
    end


    return ham
end

function make_Δtsc(Δ)
    Nx, Ny = size(Δ)
    N = Nx * Ny
    Δmat = spzeros(ComplexF64, N, N)
    for ix = 1:Nx
        for iy = 1:Ny
            i = (iy - 1) * Nx + ix
            Δmat[i, i] = Δ[ix, iy]
        end
    end
    return kron(Δmat, im * σy)
end

function make_Htsc_normal(Nx, Ny, μ, BC, h, α)
    N = Nx * Ny
    Tx_plushop = make_x_plushop(Nx, Ny, BC)
    Tx_minushop = make_x_minushop(Nx, Ny, BC)
    Ty_plushop = make_y_plushop(Nx, Ny, BC)
    Ty_minushop = make_y_minushop(Nx, Ny, BC)
    HN = kron(sparse(I, N, N) * (-μ), σ0)
    HN += kron(sparse(I, N, N) * (-h), σz) #Zeeman magnetic field
    t = 1.0

    HN += kron(-t * (Tx_plushop + Tx_minushop + Ty_plushop + Ty_minushop), σ0)

    Hax = kron((α / (2im)) * (Tx_plushop - Tx_minushop), σy)
    HN += Hax
    Hay = kron((α / (2im)) * (Ty_plushop - Ty_minushop), σx)
    HN += Hay

    return HN
end

function make_Htsc_sc(Nx, Ny, μ, Δ, BC, h, α)
    HN = make_Htsc_normal(Nx, Ny, μ, BC, h, α)
    matΔ = make_Δtsc(Δ)
    H = [
        HN matΔ
        matΔ' -conj.(HN)
    ]
    return H
end

function testTSC2()
    μ = 3.5
    Δ0 = 1#0.35


    T = 0.01
    U = -5.6
    h = 2
    α = 1

    Nxs = [16, 20, 24, 30, 32, 36, 40, 44, 48, 52, 56, 60, 64, 70, 76, 80, 82, 90, 100, 120, 128, 160, 320]

    isxPBC = false
    isyPBC = false
    num_internal_degree = 2



    fp = open("NxdepNq100.txt", "w")
    for Nx in Nxs
        Ny = Nx
        ix = Nx ÷ 2
        iy = Ny ÷ 2
        i = (iy - 1) * Nx + ix
        Δs = Δ0 * ones(Nx * Ny)
        Δsnew = similar(Δs)
        ham = make_TSC_hamiltonian(Nx, Ny, μ, Δs, h, α, isxPBC, isyPBC)
        method = "Contour"
        projector = Projector(ham; method, Nq=100)
        method = "Exact"
        projector2 = Projector(ham; method, Nq=100)

        C = real(make_C(projector, i, Nx))
        C2 = real(make_C(projector2, i, Nx))

        t1 = @elapsed real(make_C(projector, i, Nx))
        t2 = @elapsed real(make_C(projector2, i, Nx))
        println("$Nx $C $C2 $(abs(C-C2)/abs(C)) $t1 $t2")
        println(fp, "$Nx $C $C2 $(abs(C-C2)/abs(C)) $t1")

    end
    close(fp)


end


function testTSC()
    Nx = 16 * 2
    Ny = 16 * 2
    μ = 3.5
    Δ0 = 1#0.35
    Δs = Δ0 * ones(Nx * Ny)
    Δsnew = similar(Δs)
    T = 0.01
    U = -5.6
    h = 2
    α = 1

    isxPBC = false
    isyPBC = false
    BC = "OBC"
    num_internal_degree = 2
    #H = make_Htsc_sc(Nx, Ny, μ, Δ0 * ones(Nx, Ny), BC, h, α)
    #ham = Hamiltonian(H, Nx * Ny; num_internal_degree)
    ham = make_TSC_hamiltonian(Nx, Ny, μ, Δs, h, α, isxPBC, isyPBC)

    #diffH = H - ham.matrix
    #display(diffH)
    #=

    n, _ = size(H)
    for i = 1:n
        for j = 1:n
            if abs(H[i, j]) > 1e-5
                if abs(H[i, j] - ham.matrix[i, j]) > 1e-6
                    println("$i $j $(H[i,j]) $(ham.matrix[i,j])")
                end
            end
            if abs(ham.matrix[i, j]) > 1e-5
                if abs(H[i, j] - ham.matrix[i, j]) > 1e-6
                    println("$i $j $(H[i,j]) $(ham.matrix[i,j])")
                end
            end
        end
    end
    =#
    #return

    #projector = Projector(ham)
    method = "Contour"
    projector = Projector(ham; method, Nq=10)
    method = "Exact"
    #projector2 = Projector(ham; method)
    debugmode = false
    fp = open("Nqdep.txt", "w")

    ix = Nx ÷ 2
    iy = Ny ÷ 2
    i = (iy - 1) * Nx + ix
    for Nq in [10, 50, 100, 150, 200, 300, 400, 500]
        method = "Contour"
        projector = Projector(ham; method, Nq)
        C = make_C(projector, i, Nx; debugmode)
        println("$Nq $(real(C))")
        println(fp, "$Nq $(real(C))")
        @time C = make_C(projector, i, Nx; debugmode)

    end
    close(fp)

    C = make_C(projector, i, Nx; debugmode)
    println("$i $C")

    projector2 = Projector(ham; method)
    C = make_C(projector2, i, Nx; debugmode)
    println("$i $C")
    return



    fp = open("test2.txt", "w")
    for ix = 1:Nx
        for iy = 1:Ny
            i = (iy - 1) * Nx + ix
            #for i = 1:Nx*Ny
            C = make_C(projector, i, Nx; debugmode)
            println("$i $C")
            println(fp, "$ix $iy $(real(C))")
            #C = make_C(projector2, i, Nx; debugmode)
            #println("$i $C")
            #println(fp, "$ix $iy $(real(C))")
        end
        println(fp, "\t")
    end
    close(fp)

end

function test5()

    μ = -1.5
    Δ = 0.5

    Nx = 16 * 4 * 2
    Ny = 16 * 4 * 2
    N = Nx * Ny
    Δ = 0.5
    ham = Hamiltonian(N, isSC=true)
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


            jx = ix - 1
            jx += ifelse(jx < 1, Nx, 0)
            jy = iy
            j = (jy - 1) * Nx + jx
            cj = FermionOP(j)
            ham += t * (ci' * cj - ci * cj')


            jy = iy + 1
            jy += ifelse(jy > Ny, -Ny, 0)
            jx = ix
            j = (jy - 1) * Nx + jx
            cj = FermionOP(j)
            ham += t * (ci' * cj - ci * cj')


            jy = iy - 1
            jy += ifelse(jy < 1, Ny, 0)
            jx = ix
            j = (jy - 1) * Nx + jx
            cj = FermionOP(j)
            ham += t * (ci' * cj - ci * cj')


            j = i
            cj = FermionOP(j)
            ham += -μ * (ci' * cj - ci * cj')


            ham += Δ * ci' * cj' + Δ * ci * cj

        end
    end


    T = 0.1
    m = Meanfields_solver(ham, T, method="Chebyshev", nmax=200)
    m2 = Meanfields_solver(ham, T)
    c1up = FermionOP(1, 1)
    c1down = FermionOP(1, 2)
    Gij0 = calc_meanfields(m, c1up, c1down) #<c1up c1down>
    @btime calc_meanfields($m, $c1up, $c1down) #<c1up c1down>
    println("Chebyshev: ", Gij0)
    Gij0 = calc_meanfields(m2, c1up, c1down) #<c1up c1down>
    @btime calc_meanfields($m2, $c1up, $c1down) #<c1up c1down>
    println("RSCG ", Gij0)

    m3 = Meanfields_solver(ham, T, method="Chebyshev", isLK=true, nmax=200)
    Gij0 = calc_meanfields(m3, c1up, c1down) #<c1up c1down>
    @btime calc_meanfields($m3, $c1up, $c1down) #<c1up c1down>
    println("Chebyshev: ", Gij0)
end

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
    test()
    test2()
    test3()
    ##test4()
    test5()
    #testTSC2()
end
