using QuadraticHamiltonians

function make_TSC_hamiltonian(Nx, Ny, μ, Δs, h, α, isPBC)
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
                    if isPBC
                        jx += ifelse(jx > Nx, -Nx, 0)
                        jx += ifelse(jx < 1, Nx, 0)
                    end
                    jy = iy + dy
                    if isPBC
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
                            else
                                ham += -1 * (α / (2im)) * (ci' * cj - cj * ci') * σy
                            end
                        elseif dx == 0
                            if dy == 1
                                ham += (α / (2im)) * (ci' * cj - cj * ci') * σx
                            else
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

function update!(m, Δs)
    N, _ = size(m.hamiltonian)
    for i = 1:(N÷4)
        ciup = FermionOP(i, 1)
        cidown = FermionOP(i, 2)
        update_hamiltonian!(m, -Δs[i], ciup', cidown')
        update_hamiltonian!(m, Δs[i], cidown', ciup')
        update_hamiltonian!(m, Δs[i], ciup, cidown)
        update_hamiltonian!(m, -Δs[i], cidown, ciup)
    end
end



function main()
    Nx = 16
    Ny = 16
    μ = 3.5
    Δ0 = 3
    Δs = Δ0 * ones(Nx * Ny)
    Δsnew = similar(Δs)
    T = 0.01
    U = -5.6
    h = 1
    α = 1

    isPBC = false
    ham = make_TSC_hamiltonian(Nx, Ny, μ, Δs, h, α, isPBC)
    m = Meanfields_solver(ham, T)

    for it = 1:100
        @time Gs = map(i -> calc_meanfields(m, FermionOP(i, 1), FermionOP(i, 2)), 1:Nx*Ny)
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