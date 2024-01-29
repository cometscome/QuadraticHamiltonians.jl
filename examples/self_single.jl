using QuadraticHamiltonians

function make_hamiltonian(Nx, Ny, μ, Δs)
    N = Nx * Ny
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
            ham += Δs[i] * ci' * ci' + Δs[i] * ci * ci
        end
    end
    return ham
end

function update!(m, Δs)
    N, _ = size(m.hamiltonian)
    for i = 1:(N÷2)
        ci = FermionOP(i)
        update_hamiltonian!(m, Δs[i], ci', ci')
        update_hamiltonian!(m, Δs[i], ci, ci)
    end
end

function main2()

    Nx = 24
    Ny = 24
    μ = -0.2
    Δ0 = 0.5
    Δs = Δ0 * ones(Nx * Ny)
    Δsnew = similar(Δs)
    T = 0.02
    U = -2
    aa = 10
    nmax = 200
    isLK = false

    ham = make_hamiltonian(Nx, Ny, μ, Δs)
    m = Meanfields_solver(ham, T, method="Chebyshev"; aa, nmax, isLK)

    for it = 1:100
        @time Gs = map(i -> calc_meanfields(m, FermionOP(i), FermionOP(i)), 1:Nx*Ny)
        Δsnew .= real(U * Gs)
        res = sum(abs.(Δsnew .- Δs)) / sum(abs.(Δs))
        println("$(it)-th $(Δsnew[1]) $res")
        if res < 1e-4
            break
        end
        update!(m, Δsnew)
        Δs .= Δsnew
    end
    ham = get_hamiltonian(m)

end


function main()

    Nx = 24
    Ny = 24
    μ = -0.2
    Δ0 = 0.5
    Δs = Δ0 * ones(Nx * Ny)
    Δsnew = similar(Δs)
    T = 0.02
    U = -2

    ham = make_hamiltonian(Nx, Ny, μ, Δs)
    m = Meanfields_solver(ham, T)

    for it = 1:100
        @time Gs = map(i -> calc_meanfields(m, FermionOP(i), FermionOP(i)), 1:Nx*Ny)
        Δsnew .= real(U * Gs)
        res = sum(abs.(Δsnew .- Δs)) / sum(abs.(Δs))
        println("$(it)-th $(Δsnew[1]) $res")
        if res < 1e-4
            break
        end
        update!(m, Δsnew)
        Δs .= Δsnew
    end

end
main()
main2()