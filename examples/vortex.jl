using QuadraticHamiltonians
using Plots

function main()
    μ = -1.5
    Nx = 128
    Ny = 128
    N = Nx * Ny
    ξ = 2
    Δ0 = 0.2
    Δxy(x, y) = Δ0 * tanh(sqrt(x^2 + y^2) / ξ) * exp(im * atan(y, x))
    calc_x(ix) = ix - (Nx + 1) / 2
    calc_y(iy) = iy - (Ny + 1) / 2

    ham = Hamiltonian(ComplexF64, N, isSC=true)
    t = -1
    hops = [(+1, 0), (-1, 0), (0, +1), (0, -1)]
    for ix = 1:Nx
        for iy = 1:Ny
            i = (iy - 1) * Nx + ix
            ci = FermionOP(i)
            for (dx, dy) in hops
                jx = ix + dx
                jy = iy + dy
                if 1 <= jx <= Nx && 1 <= jy <= Ny
                    j = (jy - 1) * Nx + jx
                    cj = FermionOP(j)
                    ham += -1 * (ci' * cj - cj * ci')
                end
            end
            ham += -μ * (ci' * ci - ci * ci')
            x = calc_x(ix)
            y = calc_y(iy)
            ham += Δxy(x, y) * ci' * ci' + Δxy(x, y)' * ci * ci
        end
    end

    num = 200
    zs = zeros(ComplexF64, num)
    ene = range(-1, 1, length=num)
    eta = 0.05
    for i = 1:num
        zs[i] = ene[i] + im * eta
    end
    npoints = 20
    rhoz = zeros(num, npoints)
    f = plot()
    for ip = 0:npoints-1
        ix = Nx ÷ 2 + ip
        iy = ix
        i = (iy - 1) * Nx + ix
        c1 = FermionOP(i)
        println(ip)
        @time Gij = calc_greenfunction(ham, zs, c1', c1)
        rhoz[:, ip+1] = -imag.(Gij) / π
        plot!(f, rhoz[:, ip+1], label=ip)

    end
    savefig(f, "ldosv.png")
end

main()