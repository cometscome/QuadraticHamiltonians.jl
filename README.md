# QuadraticHamiltonians.jl

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://cometscome.github.io/QuadraticHamiltonians.jl/stable/)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://cometscome.github.io/QuadraticHamiltonians.jl/dev/)
[![Build Status](https://github.com/cometscome/QuadraticHamiltonians.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/cometscome/QuadraticHamiltonians.jl/actions/workflows/CI.yml?query=branch%3Amain)


This is an alpha version. 

This calculates meanfields of the Hamiltonian with quadratic terms in superconductors. 
With the use of this package, you do not have to think about a matrix of the Hamiltonian that you want to consider. To construct the Hamiltonian, just add the operators. 


We implemented two kinds of methods. 

-  Chebyshev Polynomial method: Yuki Nagai, Yukihiro Ota, and Masahiko Machida, Efficient Numerical Self-Consistent Mean-Field Approach for Fermionic Many-Body Systems by Polynomial Expansion on Spectral Density, [J. Phys. Soc. Jpn. 81, 024710 (2012)](https://journals.jps.jp/doi/10.1143/JPSJ.81.024710)
-  Reduced-Shifted Conjugate-Gradient Method method: Yuki Nagai et al., Reduced-Shifted Conjugate-Gradient Method for a Green’s Function: Efficient Numerical Approach in a Nano-Structured Superconductor, [J. Phys. Soc. Jpn. 86, 014708 (2017)](https://journals.jps.jp/doi/10.7566/JPSJ.86.014708). In the RSCG method, we use the sparse modeling technique proposed by [SparseIR.jl](https://github.com/SpM-lab/SparseIR.jl). 

And also we implemented the following method:
- 1. LK-BdG solver: Yuki Nagai, N-independent Localized Krylov–Bogoliubov-de Gennes Method: Ultra-fast Numerical Approach to Large-scale Inhomogeneous Superconductors [J. Phys. Soc. Jpn. 89, 074703 (2020) ](https://journals.jps.jp/doi/10.7566/JPSJ.89.074703)
 
Please note that this is an alpha version. 

We consider the Hamiltonian only with quadratic form terms in real space.
```math
\hat{H} \equiv \sum_{i j} \sum_{\alpha \beta} H_{ij} c_{i \alpha}^{\dagger} c_{j \beta}
```
or
```math
\hat{H}_{\rm BdG} \equiv \sum_{i j} \sum_{\alpha \beta} H_{ij} c_{i \alpha}^{\dagger} c_{j \beta} + \sum_{i j} \sum_{\alpha \beta} \bar{H}_{ij} c_{i \alpha} c_{j \beta}^{\dagger} + + \sum_{i j} \sum_{\alpha \beta} \Delta_{ij} c_{i \alpha}^{\dagger} c_{j \beta}^{\dagger} +  \sum_{i j} \sum_{\alpha \beta} \bar{\Delta}_{ij} c_{i \alpha} c_{j \beta} 
```
if the keyword ```isSC``` is ```true```. Here $i$ and $j$ are indices of real-space grid and $\alpha$ and $\beta$ are indices of internal degree of freedom (like spins,orbitals or bands). 
If you want to consider a superconducting BdG Hamiltonian, you need the relations:
```math
\begin{align}
\bar{H}_{ij} &= - H_{ij}^{\ast}, \\
\bar{\Delta}_{ij} &= \Delta_{ji}^{\ast}.
\end{align}
``` 

Now this package can 
- construct the matrix of the Hamiltonian.
- display the creation and anihilation operators in the Hamiltonian.
-  calculate meanfields defined as $\langle c_i^{\dagger} c_j \rangle$, $\langle c_i^{\dagger} c_j^{\dagger} \rangle$ , $\langle c_i c_j^{\dagger} \rangle$ and $\langle c_i c_j \rangle$. You can choose the Chebyshev polynomial method or RSCG method for calculating meanfields.  Now only a Hermitian Hamiltonian is supported. 
-  calculate Green's functions defined as $G_{ij}(z) \equiv \langle c_i^{\dagger} [z \hat{I} - \hat{H}]^{-1} c_j \rangle = \left[ [z \hat{I} - \hat{H}]^{-1} \right]_{ij}$. Here $z$ is a complex frequency. Now only a Hermitian Hamiltonian is supported. 


# Install
```
add https://github.com/cometscome/QuadraticHamiltonians.jl
```
# How to construct Hamiltonians
## Normal states with single orbital
Let us consider a normal state 4-site Hamiltonian
```
julia> H = Hamiltonian(4)
---------------------------------
Hamiltonian: 
Num. of sites: 4
Num. of internal degree of freedom: 1
H = 
---------------------------------
```
We define operators
```julia
julia> c1 = FermionOP(1)
+C_{1,1}

julia> c2 = FermionOP(2)
+C_{2,1}
```

Then we can add operators
```julia
julia> H += 2.0*c1'*c1 + 3*c2'*c2 + 2*c1'*c2 + 2*(c1'*c2)'
---------------------------------
Hamiltonian: 
Num. of sites: 4
Num. of internal degree of freedom: 1
H = +2.0C_{1,1}^+C_{1,1} +2.0C_{2,1}^+C_{1,1} +2.0C_{1,1}^+C_{2,1} +3.0C_{2,1}^+C_{2,1} 
---------------------------------
```
We can regard H as a matrix. So you can do like 
```julia
julia> x = rand(4)
4-element Vector{Float64}:
 0.18622544459032153
 0.36202722147302
 0.9038045102532898
 0.8497382867289013

julia> H*x
4-element Vector{Float64}:
 1.096505332126683
 1.4585325535997031
 0.0
 0.0
```

If we want to consider the Hamiltonian whose elements are complex values, we define 
```julia
julia> H = Hamiltonian(ComplexF64,4)
---------------------------------
Hamiltonian: 
Num. of sites: 4
Num. of internal degree of freedom: 1
H = 
---------------------------------
```
And we can add operators:
```julia
julia> H += 2.0*c1'*c1 + 3*c2'*c2 + exp(0.2*im)*c1'*c2 + (exp(0.2*im)*c1'*c2)' 
---------------------------------
Hamiltonian: 
Num. of sites: 4
Num. of internal degree of freedom: 1
H = +2.0C_{1,1}^+C_{1,1} +(0.9800665778412416 - 0.19866933079506122im)C_{2,1}^+C_{1,1} +(0.9800665778412416 + 0.19866933079506122im)C_{1,1}^+C_{2,1} +3.0C_{2,1}^+C_{2,1} 
---------------------------------
```

## Normal states with multi orbitals
If we want to consider internal degree of freedom such as spins, orbitals and bands, we can define 
```julia
julia> H = Hamiltonian(4,num_internal_degree=2)
---------------------------------
Hamiltonian: 
Num. of sites: 4
Num. of internal degree of freedom: 2
H = 
---------------------------------
```
Operators are defined as 
```julia
julia> c1up = FermionOP(1,1)
+C_{1,1}

julia> c1down = FermionOP(1,2)
+C_{1,2}

julia> c2down = FermionOP(2,2)
+C_{2,2}
```

## Superconducting states
Let us use ```isSC``` to consider superconducting states. 
```julia
julia> H = Hamiltonian(4,isSC=true)
---------------------------------
Hamiltonian: 
Num. of sites: 4
Num. of internal degree of freedom: 1
Superconducting state
H = 
---------------------------------
```

```julia
julia> H += -1*(c1'*c1 - c1*c1') + 0.2*c2'c2' + (0.2*c2'*c2')'
---------------------------------
Hamiltonian: 
Num. of sites: 4
Num. of internal degree of freedom: 1
Superconducting state
H = -1.0C_{1,1}^+C_{1,1} +0.2C_{2,1}C_{2,1} +C_{1,1}C_{1,1}^+ +0.2C_{2,1}^+C_{2,1}^+ 
---------------------------------
```
We can have a matrix form. 
```julia
julia> H_sp = construct_matrix(H)
8×8 SparseArrays.SparseMatrixCSC{Float64, Int64} with 4 stored entries:
 -1.0   ⋅    ⋅    ⋅    ⋅    ⋅    ⋅    ⋅ 
   ⋅    ⋅    ⋅    ⋅    ⋅   0.2   ⋅    ⋅ 
   ⋅    ⋅    ⋅    ⋅    ⋅    ⋅    ⋅    ⋅ 
   ⋅    ⋅    ⋅    ⋅    ⋅    ⋅    ⋅    ⋅ 
   ⋅    ⋅    ⋅    ⋅   1.0   ⋅    ⋅    ⋅ 
   ⋅   0.2   ⋅    ⋅    ⋅    ⋅    ⋅    ⋅ 
   ⋅    ⋅    ⋅    ⋅    ⋅    ⋅    ⋅    ⋅ 
   ⋅    ⋅    ⋅    ⋅    ⋅    ⋅    ⋅    ⋅ 
```
## Tight-binding model example

```julia
using QuadraticHamiltonians
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
test2()

```

The output is 

```
---------------------------------
Hamiltonian: 
Num. of sites: 16
Num. of internal degree of freedom: 1
H = -1.0C_{1,1}^+C_{2,1} -1.0C_{1,1}^+C_{16,1} -1.0C_{2,1}^+C_{3,1} -1.0C_{2,1}^+C_{1,1} -1.0C_{3,1}^+C_{4,1} -1.0C_{3,1}^+C_{2,1} -1.0C_{4,1}^+C_{5,1} -1.0C_{4,1}^+C_{3,1} -1.0C_{5,1}^+C_{6,1} -1.0C_{5,1}^+C_{4,1} -1.0C_{6,1}^+C_{7,1} -1.0C_{6,1}^+C_{5,1} -1.0C_{7,1}^+C_{8,1} -1.0C_{7,1}^+C_{6,1} -1.0C_{8,1}^+C_{9,1} -1.0C_{8,1}^+C_{7,1} -1.0C_{9,1}^+C_{10,1} -1.0C_{9,1}^+C_{8,1} -1.0C_{10,1}^+C_{11,1} -1.0C_{10,1}^+C_{9,1} -1.0C_{11,1}^+C_{12,1} -1.0C_{11,1}^+C_{10,1} -1.0C_{12,1}^+C_{13,1} -1.0C_{12,1}^+C_{11,1} -1.0C_{13,1}^+C_{14,1} -1.0C_{13,1}^+C_{12,1} -1.0C_{14,1}^+C_{15,1} -1.0C_{14,1}^+C_{13,1} -1.0C_{15,1}^+C_{16,1} -1.0C_{15,1}^+C_{14,1} -1.0C_{16,1}^+C_{1,1} -1.0C_{16,1}^+C_{15,1} 
---------------------------------
16×16 SparseArrays.SparseMatrixCSC{Float64, Int64} with 32 stored entries:
⎡⠪⡢⡀⠀⠀⠀⠀⠈⎤
⎢⠀⠈⠪⡢⡀⠀⠀⠀⎥
⎢⠀⠀⠀⠈⠪⡢⡀⠀⎥
⎣⡀⠀⠀⠀⠀⠈⠪⡢⎦
```

## 2-dimensional s-wave superconductor
```julia
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
```

# How to calculate meanfields
You can calculate meanfields. 
To calculate meanfields, we define the solver: 
```julia
T = 0.1
m = Meanfields_solver(ham, T)
```
Then the RSCG method will be used. If we want to use the Chebyshev polynomial method, we define 
```julia
m = Meanfields_solver(ham, T, method="Chebyshev", nmax=200)
```
If we want to use the LK-Chebyshev polynomial method, we define 
```julia
m = Meanfields_solver(ham, T, method="Chebyshev", isLK=true, nmax=200)
```

For example, meanfield $\langle c_1 c_1 \rangle$ is calculated by
```julia
c1 = FermionOP(1)
Gij0 = calc_meanfields(m, c1, c1) #<c1 c1>
```

If we want to update the Hamiltonian, we can use 
```julia
update_hamiltonian!(m, value, ci', ci')
update_hamiltonian!(m, value, ci, ci)
```

If we want to get the Hamiltonian, we can use 
```julia
ham = get_hamiltonian(m)
```

# How to calculate Green's functions or local density of states


# Examples

## s-wave superconductor

```julia
using QuadraticHamiltonians
using Plots
using BenchmarkTools

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
```

The output is 
```
The solver is the Chebyshev polynomial method
  27.004 ms (10 allocations: 1.00 MiB)
Chebyshev: -0.1761175169590693
The solver is the RSCG
num. of Matsubara freqs. 24
  33.507 ms (6412 allocations: 37.76 MiB)
RSCG -0.17610594510770422 - 4.9649942910811176e-15im
The solver is the Chebyshev polynomial method
  20.113 ms (42 allocations: 2.13 MiB)
Chebyshev with LKvectors: -0.17611751734430567
```

![ldos](https://github.com/cometscome/QuadraticHamiltonians.jl/assets/21115243/bde82aa4-0506-4d9e-8742-19dd61706b76)

## self-consistent gap calculation
```julia
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

function update_hamiltonian!(ham, Δold, Δs)
    N, _ = size(ham)
    for i = 1:(N÷2)
        ci = FermionOP(i)
        ham += (Δs[i] - Δold[i]) * ci' * ci' + (Δs[i] - Δold[i]) * ci * ci
    end
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
        update_hamiltonian!(ham, Δs, Δsnew)
        Δs .= Δsnew
    end

end
main()
```

The output is 
```
The solver is the RSCG
num. of Matsubara freqs. 34
  3.081818 seconds (14.92 M allocations: 3.000 GiB, 7.47% gc time, 28.05% compilation time)
1-th 0.43401688579960984 0.13196622709719913
  2.425116 seconds (11.58 M allocations: 3.088 GiB, 4.43% gc time)
2-th 0.4020473310731843 0.07365970445490454
  2.556391 seconds (11.96 M allocations: 3.262 GiB, 4.12% gc time)
3-th 0.3853252114191848 0.04159241583023033
  2.628881 seconds (12.24 M allocations: 3.358 GiB, 4.02% gc time)
4-th 0.3762161124243843 0.02364002840649464
  2.708366 seconds (12.36 M allocations: 3.398 GiB, 4.26% gc time)
5-th 0.37114204479276297 0.013487108372129082
  2.765002 seconds (12.34 M allocations: 3.397 GiB, 4.81% gc time)
6-th 0.36828006208588193 0.0077112866961116305
  2.743782 seconds (12.56 M allocations: 3.476 GiB, 4.05% gc time)
7-th 0.36665432701604384 0.0044143990509399635
  2.723989 seconds (12.54 M allocations: 3.490 GiB, 3.88% gc time)
8-th 0.36572710332447894 0.0025288823568850502
  2.744749 seconds (12.73 M allocations: 3.520 GiB, 4.05% gc time)
9-th 0.3651970461434803 0.001449321833805112
  2.770678 seconds (12.65 M allocations: 3.515 GiB, 4.28% gc time)
10-th 0.3648936376188999 0.0008308059042921653
  2.741237 seconds (12.65 M allocations: 3.515 GiB, 4.07% gc time)
11-th 0.36471983333293023 0.00047631387916426306
  2.771782 seconds (12.67 M allocations: 3.514 GiB, 4.30% gc time)
12-th 0.36462022832989466 0.0002730987087264538
  2.748457 seconds (12.67 M allocations: 3.515 GiB, 4.08% gc time)
13-th 0.36456313174463706 0.00015659133320500504
  2.820919 seconds (12.69 M allocations: 3.518 GiB, 4.51% gc time)
14-th 0.36453039805715926 8.978993383823577e-5
```

## s-wave superconductor with spins

```julia
using QuadraticHamiltonians
using Plots
using BenchmarkTools
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
main2()
```

The output is 
```
The solver is the Chebyshev polynomial method
  49.120 ms (10 allocations: 2.00 MiB)
Chebyshev: 0.1761175169590693
The solver is the RSCG
num. of Matsubara freqs. 24
  63.067 ms (6412 allocations: 75.01 MiB)
RSCG 0.17610594510770403 + 4.955324870183428e-15im
The solver is the Chebyshev polynomial method
  21.877 ms (42 allocations: 4.25 MiB)
Chebyshev with LKvectors: 0.1761175173443056
```

![ldos_spin](https://github.com/cometscome/QuadraticHamiltonians.jl/assets/21115243/bd19da9d-ab50-42bc-99df-bf0bed5d888d)
