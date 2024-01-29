using QuadraticHamiltonians
using Plots
using BenchmarkTools

function test()
    H = Hamiltonian(16)
    display(H)
    c1 = FermionOP(1)
    c3 = FermionOP(3)
    c4 = FermionOP(4)
    display(0.3 * c4' * c3)
    display(im * c1' * c4)
    H += -1 * c1' * c1 + 2 * c3' * c1 + 0.3 * c4' * c3
    display(H)
    H -= 2 * c3' * c1 + 0.3 * c4' * c3
    display(H)
end
test()

function main()
    μ = -1.5
    Nx = 128
    Ny = 128
    N = Nx * Ny
    Δ = 0.5
    ham = Hamiltonian(N, isSC=true)
    t = -1
    hops = [(+1, 0), (-1, 0), (0, +1), (0, -1)]
    for ix = 1:Nx
        for iy = 1:Ny
            i = (iy - 1) * Nx + ix
            ci = FermionOP(i)
            for (dx, dy) in hops
                jx = ix + dx
                jx += ifelse(jx > Nx, -Nx, 0)
                jx += ifelse(jx < 1, Nx, 0)
                jy = iy + dy
                jy += ifelse(jy > Ny, -Ny, 0)
                jy += ifelse(jy < 1, Ny, 0)
                j = (jy - 1) * Nx + jx
                cj = FermionOP(j)
                ham += -1 * (ci' * cj - ci * cj')
            end
            ham += -μ * (ci' * ci - ci * ci')
            ham += Δ * ci' * ci' + Δ * ci * ci
        end
    end

    num = 200
    zs = zeros(ComplexF64, num)
    ene = range(-6, 6, length=num)
    eta = 0.05
    for i = 1:num
        zs[i] = ene[i] + im * eta
    end
    ix = Nx ÷ 2
    iy = ix
    i = (iy - 1) * Nx + ix
    c1 = FermionOP(i)
    Gij = calc_greenfunction(ham, zs, c1', c1)
    plot(ene, -imag.(Gij) / π)
    savefig("ldos.png")

    T = 0.1
    m = Meanfields_solver(ham, T, method="Chebyshev", nmax=200)
    c1up = FermionOP(1)
    c1down = FermionOP(1)
    Gij0 = calc_meanfields(m, c1up, c1down) #<c1up c1down>
    @btime calc_meanfields($m, $c1up, $c1down) #<c1up c1down>
    println("Chebyshev: ", Gij0)
    m2 = Meanfields_solver(ham, T)
    Gij0 = calc_meanfields(m2, c1up, c1down) #<c1up c1down>
    @btime calc_meanfields($m2, $c1up, $c1down) #<c1up c1down>
    println("RSCG ", Gij0)

    m3 = Meanfields_solver(ham, T, method="Chebyshev", isLK=true, nmax=200)
    Gij0 = calc_meanfields(m3, c1up, c1down) #<c1up c1down>
    @btime calc_meanfields($m3, $c1up, $c1down) #<c1up c1down>
    println("Chebyshev with LKvectors: ", Gij0)

end


function main2()
    μ = -1.5
    Nx = 128
    Ny = 128
    N = Nx * Ny
    Δ = 0.5
    ham = Hamiltonian(N, num_internal_degree=2, isSC=true)
    t = -1
    hops = [(+1, 0), (-1, 0), (0, +1), (0, -1)]
    for ix = 1:Nx
        for iy = 1:Ny
            i = (iy - 1) * Nx + ix
            for ispin = 1:2
                ci = FermionOP(i, ispin)
                jspin = ispin
                for (dx, dy) in hops
                    jx = ix + dx
                    jx += ifelse(jx > Nx, -Nx, 0)
                    jx += ifelse(jx < 1, Nx, 0)
                    jy = iy + dy
                    jy += ifelse(jy > Ny, -Ny, 0)
                    jy += ifelse(jy < 1, Ny, 0)
                    j = (jy - 1) * Nx + jx
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
    end

    num = 200
    zs = zeros(ComplexF64, num)
    ene = range(-6, 6, length=num)
    eta = 0.05
    for i = 1:num
        zs[i] = ene[i] + im * eta
    end
    ix = Nx ÷ 2
    iy = ix
    i = (iy - 1) * Nx + ix
    c1 = FermionOP(i, 1)
    Gij = calc_greenfunction(ham, zs, c1', c1)
    plot(ene, -imag.(Gij) / π)
    savefig("ldos_spin.png")

    T = 0.1
    m = Meanfields_solver(ham, T, method="Chebyshev", nmax=200)
    c1up = FermionOP(1, 1)
    c1down = FermionOP(1, 2)
    Gij0 = calc_meanfields(m, c1up, c1down) #<c1up c1down>
    @btime calc_meanfields($m, $c1up, $c1down) #<c1up c1down>
    println("Chebyshev: ", Gij0)

    m2 = Meanfields_solver(ham, T)
    Gij0 = calc_meanfields(m2, c1up, c1down) #<c1up c1down>
    @btime calc_meanfields($m2, $c1up, $c1down) #<c1up c1down>
    println("RSCG ", Gij0)

    m3 = Meanfields_solver(ham, T, method="Chebyshev", isLK=true, nmax=200)
    Gij0 = calc_meanfields(m3, c1up, c1down) #<c1up c1down>
    @btime calc_meanfields($m3, $c1up, $c1down) #<c1up c1down>
    println("Chebyshev with LKvectors: ", Gij0)


end
#main()
#main2()