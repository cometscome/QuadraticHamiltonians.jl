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
-  calculate Green's functions defined as $G_{ij}(z) \equiv \langle c_i^{\dagger} [z \hat{I} - \hat{H}]^{-1} c_j \rangle = \left[ [z \hat{I} - \hat{H}]^{-1} \right]_{ij}$. Here $z$ is a complex frequency. Now only a Hermitian Hamiltonian is supported. The local density of states can be caluculated by the imaginary part of diagonal elements of the Green's function. 


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
There are two different ways for setting a value: 
```julia
H[2, 1] = 21
H[c3', c1] = 42
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

To get elements of the Hamiltonian, you can do like
```julia
a = H[1,2]
```
or 
```julia
c1 = FermionOP(1)
c2 = FermionOP(2)
a = H[c1',c2]
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
            ham += -1 * (ci' * cj - cj * ci')
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
If the Hamiltonian is expressed as the sum of the quadratic terms, the Green's function defined as 
```math
\begin{align}
G_{ij}(z) \equiv \langle c_i^{\dagger} [z \hat{I} - \hat{H}]^{-1} c_j \rangle
\end{align}
```
can be calculated by 
```math
\begin{align}
G_{ij}(z) = \left[ (z \hat{I} - \hat{H})^{-1} \right]_{ij}
\end{align}
```
The local density of states at site $i$ can be calculated by the Green's function: 
```math
\begin{align}
N_i(\omega) = -\frac{1}{\pi} {\rm Im} \: \lim_{\eta \rightarrow 0+} G_{ii}(\omega + i \eta) 
\end{align}
```

In this package, the Green's function can be calculated by 
```julia
z = 0.1 + 0.2im
i = 1
c1 = FermionOP(i)
Gij = calc_greenfunction(ham, z, c1', c1)
```
with the use of the RSCG method. 
Since the RSCG method can consider many complex frequencies simultaneously, we can calculate the Green's function with different frequencies $z_k$: 
```julia
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
```



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
                ham += -1 * (ci' * cj - cj * ci')
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
                ham += -1 * (ci' * cj - cj * ci')
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
```

The output is 
```
The solver is the RSCG
num. of Matsubara freqs. 34
  1.357527 seconds (8.97 M allocations: 1.225 GiB, 7.59% gc time, 46.78% compilation time)
1-th 0.43392863318462976 0.13214274366572432
  0.827926 seconds (5.99 M allocations: 1.104 GiB, 5.64% gc time)
2-th 0.40181627897040467 0.07400376230129402
  0.881787 seconds (6.16 M allocations: 1.144 GiB, 5.09% gc time)
3-th 0.38493556921008637 0.04201103156997562
  0.884848 seconds (6.32 M allocations: 1.180 GiB, 5.40% gc time)
4-th 0.37568124638585454 0.024041167724526982
  0.891520 seconds (6.36 M allocations: 1.192 GiB, 5.27% gc time)
5-th 0.3704883477225273 0.013822612441196765
  0.889258 seconds (6.38 M allocations: 1.215 GiB, 4.45% gc time)
6-th 0.3675358557663468 0.007969231662928523
  0.904341 seconds (6.44 M allocations: 1.212 GiB, 5.19% gc time)
7-th 0.36584446964213757 0.004601915465106033
  0.946217 seconds (6.52 M allocations: 1.231 GiB, 5.29% gc time)
8-th 0.36487138389829876 0.0026599130537329862
  0.903315 seconds (6.49 M allocations: 1.224 GiB, 4.95% gc time)
9-th 0.36431012784136424 0.001538238218849528
  0.905315 seconds (6.52 M allocations: 1.232 GiB, 4.69% gc time)
10-th 0.3639859317169275 0.0008898496406517458
  0.904538 seconds (6.52 M allocations: 1.232 GiB, 4.80% gc time)
11-th 0.3637985531159851 0.0005148569471265923
  0.917017 seconds (6.56 M allocations: 1.240 GiB, 4.78% gc time)
12-th 0.36369015468817556 0.0002979221928351472
  0.945926 seconds (6.54 M allocations: 1.236 GiB, 5.33% gc time)
13-th 0.36362744231704996 0.00017239402356417755
  0.922917 seconds (6.50 M allocations: 1.225 GiB, 5.31% gc time)
14-th 0.3635911839621331 9.976804879798009e-5
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
                    ham += -1 * (ci' * cj - cj * ci')
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

## Topological s-wave superconductor with the Zeeman term and the Rashba spin-orbit term

```julia
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
```

The output is 
```
The solver is the RSCG
num. of Matsubara freqs. 40
  1.079232 seconds (6.06 M allocations: 707.212 MiB, 16.58% gc time, 64.22% compilation time)
1-th -1.873023776160689 1.6271730164203106
  0.601234 seconds (4.04 M allocations: 846.823 MiB, 4.71% gc time)
2-th -1.4418945612222904 0.20965023882070202
  0.918474 seconds (5.74 M allocations: 1.242 GiB, 4.02% gc time)
3-th -1.2113675935614676 0.13563888454719064
  1.213072 seconds (7.08 M allocations: 1.605 GiB, 4.38% gc time)
4-th -1.0683848898677126 0.09725266043470575
  1.519615 seconds (8.41 M allocations: 1.988 GiB, 4.20% gc time)
5-th -0.9719173169890015 0.07520905143665893
  1.774157 seconds (9.26 M allocations: 2.261 GiB, 4.52% gc time)
6-th -0.9031937676953646 0.0614030514203067
  1.960240 seconds (9.93 M allocations: 2.480 GiB, 4.87% gc time)
7-th -0.8522993398753158 0.052088904542545035
  2.113638 seconds (10.48 M allocations: 2.672 GiB, 4.82% gc time)
8-th -0.8134640044396014 0.04536935661275386
  2.238472 seconds (10.90 M allocations: 2.814 GiB, 5.02% gc time)
9-th -0.7830935035525877 0.040207500802570885
  2.127559 seconds (10.55 M allocations: 2.707 GiB, 4.57% gc time)
10-th -0.7588369948792518 0.03600999926600642
  2.334611 seconds (11.38 M allocations: 2.977 GiB, 4.48% gc time)
11-th -0.7391008228944542 0.032434049930203365
  2.442818 seconds (11.64 M allocations: 3.064 GiB, 4.63% gc time)
12-th -0.7227749038266035 0.02928518035552862
  2.496474 seconds (11.85 M allocations: 3.135 GiB, 4.41% gc time)
13-th -0.709069155872567 0.026456030788754803
  2.483085 seconds (11.97 M allocations: 3.181 GiB, 3.93% gc time)
14-th -0.6974110355699145 0.02388799764731334
  2.495291 seconds (12.05 M allocations: 3.216 GiB, 3.89% gc time)
15-th -0.6873791072315169 0.0215479761914521
  2.542803 seconds (12.16 M allocations: 3.257 GiB, 4.01% gc time)
16-th -0.6786586888786034 0.019415153943493132
  2.587446 seconds (12.25 M allocations: 3.295 GiB, 4.06% gc time)
17-th -0.671011679282689 0.017474092388203903
  2.626588 seconds (12.35 M allocations: 3.335 GiB, 4.09% gc time)
18-th -0.6642555683053231 0.01571143747168237
  2.681020 seconds (12.46 M allocations: 3.373 GiB, 4.22% gc time)
19-th -0.6582486676185052 0.014114520535718768
  2.712846 seconds (12.44 M allocations: 3.368 GiB, 4.94% gc time)
20-th -0.6528795616656254 0.012670927816306929
  2.687108 seconds (12.47 M allocations: 3.385 GiB, 4.42% gc time)
21-th -0.648059496125448 0.01136845977752193
  2.698632 seconds (12.62 M allocations: 3.432 GiB, 4.16% gc time)
22-th -0.643716834307313 0.010195258886917814
  2.774037 seconds (12.60 M allocations: 3.431 GiB, 4.52% gc time)
23-th -0.6397929998696529 0.00913994023009275
  2.804256 seconds (12.72 M allocations: 3.478 GiB, 4.61% gc time)
24-th -0.6362394570260124 0.00819171955905381
  2.775846 seconds (12.89 M allocations: 3.538 GiB, 4.08% gc time)
25-th -0.63301548265852 0.007340492709279601
  2.812354 seconds (12.98 M allocations: 3.568 GiB, 4.28% gc time)
26-th -0.6300864856968793 0.006576883266048058
  2.825387 seconds (13.03 M allocations: 3.593 GiB, 3.97% gc time)
27-th -0.6274227436671878 0.005892252532124953
  2.830010 seconds (13.08 M allocations: 3.613 GiB, 3.93% gc time)
28-th -0.6249984477604995 0.005278692772073217
  2.877537 seconds (13.10 M allocations: 3.617 GiB, 4.20% gc time)
29-th -0.6227909814572898 0.004729001813477821
  2.814024 seconds (13.14 M allocations: 3.633 GiB, 3.63% gc time)
30-th -0.6207803405763253 0.004236642241088458
  2.844506 seconds (13.15 M allocations: 3.640 GiB, 3.96% gc time)
31-th -0.6189487129414846 0.0037957051360124226
  3.005532 seconds (13.16 M allocations: 3.642 GiB, 4.95% gc time)
32-th -0.6172801374355674 0.003400858384159586
  2.982854 seconds (13.17 M allocations: 3.646 GiB, 4.76% gc time)
33-th -0.6157602345110169 0.0030473007844825273
  2.929872 seconds (13.14 M allocations: 3.638 GiB, 4.32% gc time)
34-th -0.6143759966380645 0.0027307205386831934
  2.833398 seconds (13.16 M allocations: 3.644 GiB, 3.74% gc time)
35-th -0.6131156094137203 0.0024472404597514195
  2.948779 seconds (13.18 M allocations: 3.653 GiB, 4.51% gc time)
36-th -0.6119683239104643 0.002193392792552506
  2.935843 seconds (13.21 M allocations: 3.664 GiB, 4.45% gc time)
37-th -0.6109243302390619 0.001966063888845468
  2.944080 seconds (13.22 M allocations: 3.668 GiB, 4.50% gc time)
38-th -0.6099746645043274 0.0017624669009970566
  3.032381 seconds (13.22 M allocations: 3.668 GiB, 4.98% gc time)
39-th -0.6091111276751108 0.001580108345611027
  2.912560 seconds (13.27 M allocations: 3.687 GiB, 3.91% gc time)
40-th -0.6083262105712027 0.001416755095863366
  3.070945 seconds (13.28 M allocations: 3.694 GiB, 4.85% gc time)
41-th -0.6076130419061138 0.0012704129923773262
  2.959192 seconds (13.28 M allocations: 3.693 GiB, 4.50% gc time)
42-th -0.6069653199286934 0.00113929572631485
  3.039093 seconds (13.31 M allocations: 3.700 GiB, 4.66% gc time)
43-th -0.6063772737453785 0.0010218056690919832
  3.049923 seconds (13.32 M allocations: 3.704 GiB, 4.75% gc time)
44-th -0.6058436172654165 0.0009165157705025486
  3.028901 seconds (13.31 M allocations: 3.704 GiB, 4.65% gc time)
45-th -0.6053595096728797 0.0008221475370220531
  3.014559 seconds (13.32 M allocations: 3.707 GiB, 4.55% gc time)
46-th -0.6049205193354682 0.0007375584043463334
  3.013027 seconds (13.32 M allocations: 3.709 GiB, 4.51% gc time)
47-th -0.604522590849413 0.0006617273730522861
  3.065978 seconds (13.34 M allocations: 3.713 GiB, 4.85% gc time)
48-th -0.6041620163787921 0.000593740534829054
  2.998548 seconds (13.33 M allocations: 3.712 GiB, 4.67% gc time)
49-th -0.6038354056749607 0.0005327789447378642
  3.008589 seconds (13.34 M allocations: 3.715 GiB, 4.53% gc time)
50-th -0.6035396613195223 0.00047811220217767336
  2.984136 seconds (13.34 M allocations: 3.717 GiB, 4.27% gc time)
51-th -0.603271956206971 0.0004290844742470022
  2.986102 seconds (13.34 M allocations: 3.717 GiB, 4.31% gc time)
52-th -0.6030297098728377 0.0003851109026020435
  3.008851 seconds (13.35 M allocations: 3.718 GiB, 4.27% gc time)
53-th -0.6028105687238502 0.0003456655189107658
  2.926974 seconds (13.35 M allocations: 3.720 GiB, 3.93% gc time)
54-th -0.6026123890946519 0.0003102796361298741
  3.013989 seconds (13.36 M allocations: 3.721 GiB, 4.47% gc time)
55-th -0.6024332167560249 0.00027853279022410063
  2.925957 seconds (13.35 M allocations: 3.719 GiB, 4.07% gc time)
56-th -0.6022712741026642 0.00025004790660750666
  2.990428 seconds (13.35 M allocations: 3.722 GiB, 4.25% gc time)
57-th -0.6021249422723195 0.00022448825972236396
  2.993714 seconds (13.35 M allocations: 3.719 GiB, 4.50% gc time)
58-th -0.6019927502404482 0.00020155103761652473
  2.972320 seconds (13.35 M allocations: 3.719 GiB, 4.31% gc time)
59-th -0.6018733605571956 0.0001809666539966145
  2.998017 seconds (13.36 M allocations: 3.723 GiB, 4.52% gc time)
60-th -0.6017655584431989 0.00016249232065718839
  2.958357 seconds (13.35 M allocations: 3.718 GiB, 4.23% gc time)
61-th -0.6016682408866258 0.00014591021496670423
  2.954279 seconds (13.36 M allocations: 3.721 GiB, 4.15% gc time)
62-th -0.6015804069224595 0.0001310259234340552
  2.987971 seconds (13.36 M allocations: 3.723 GiB, 4.38% gc time)
63-th -0.6015011472832793 0.00011766420738475466
  2.967308 seconds (13.35 M allocations: 3.722 GiB, 4.31% gc time)
64-th -0.6014296388255697 0.00010566954846444353
  2.952684 seconds (13.35 M allocations: 3.721 GiB, 4.15% gc time)
65-th -0.6013651365577349 9.490118655278245e-5
```

## Parallel computing
To use many cores, just use Distributed package. 

```julia
using Distributed
@everywhere using QuadraticHamiltonians
```
In the case of the topological s-wave superconductor, the code is written as 
```julia
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
        @time Gs = pmap(i -> calc_meanfields(m, FermionOP(i, 1), FermionOP(i, 2)), 1:Nx*Ny)
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
```
The only difference is ```pmap```. 

The output is 
```
julia -p 4 --project=../ self.jl
The solver is the RSCG
num. of Matsubara freqs. 40
  1.603179 seconds (1.33 M allocations: 126.028 MiB, 0.79% gc time, 25.78% compilation time)
1-th -1.873023776160689 1.6271730164203106
  0.192392 seconds (42.37 k allocations: 38.468 MiB, 2.00% gc time)
2-th -1.4418945612222904 0.20965023882070202
  0.270357 seconds (42.38 k allocations: 38.491 MiB, 0.99% gc time)
3-th -1.2113675935614676 0.13563888454719064
  0.344826 seconds (42.37 k allocations: 38.471 MiB, 0.30% gc time)
4-th -1.0683848898677126 0.09725266043470575
  0.424078 seconds (42.38 k allocations: 38.480 MiB, 0.29% gc time)
5-th -0.9719173169890015 0.07520905143665893
  0.482996 seconds (42.39 k allocations: 38.503 MiB, 0.26% gc time)
6-th -0.9031937676953646 0.0614030514203067
  0.524416 seconds (42.37 k allocations: 38.472 MiB, 0.27% gc time)
7-th -0.8522993398753158 0.052088904542545035
  0.567269 seconds (42.38 k allocations: 38.495 MiB, 0.22% gc time)
8-th -0.8134640044396014 0.04536935661275386
  0.602301 seconds (42.38 k allocations: 38.482 MiB)
9-th -0.7830935035525877 0.040207500802570885
  0.576101 seconds (42.38 k allocations: 38.477 MiB, 0.19% gc time)
10-th -0.7588369948792518 0.03600999926600642
  0.635709 seconds (42.38 k allocations: 38.475 MiB, 0.20% gc time)
11-th -0.7391008228944542 0.032434049930203365
  0.659751 seconds (42.38 k allocations: 38.492 MiB, 0.15% gc time)
12-th -0.7227749038266035 0.02928518035552862
  0.669727 seconds (42.38 k allocations: 38.491 MiB, 0.14% gc time)
13-th -0.709069155872567 0.026456030788754803
  0.679683 seconds (42.38 k allocations: 38.487 MiB, 0.14% gc time)
14-th -0.6974110355699145 0.02388799764731334
  0.687933 seconds (42.37 k allocations: 38.461 MiB, 0.14% gc time)
15-th -0.6873791072315169 0.0215479761914521
  0.697942 seconds (42.38 k allocations: 38.479 MiB, 0.14% gc time)
16-th -0.6786586888786034 0.019415153943493132
  0.723906 seconds (42.38 k allocations: 38.476 MiB)
17-th -0.671011679282689 0.017474092388203903
  0.726147 seconds (42.37 k allocations: 38.471 MiB, 0.16% gc time)
18-th -0.6642555683053231 0.01571143747168237
  0.731059 seconds (42.38 k allocations: 38.483 MiB, 0.16% gc time)
19-th -0.6582486676185052 0.014114520535718768
  0.726339 seconds (42.38 k allocations: 38.487 MiB, 0.16% gc time)
20-th -0.6528795616656254 0.012670927816306929
  0.730802 seconds (42.38 k allocations: 38.482 MiB, 0.15% gc time)
21-th -0.648059496125448 0.01136845977752193
  0.744700 seconds (42.37 k allocations: 38.462 MiB, 0.15% gc time)
22-th -0.643716834307313 0.010195258886917814
  0.739585 seconds (42.37 k allocations: 38.477 MiB, 0.14% gc time)
23-th -0.6397929998696529 0.00913994023009275
  0.749548 seconds (42.38 k allocations: 38.476 MiB, 0.16% gc time)
24-th -0.6362394570260124 0.00819171955905381
  0.756281 seconds (42.37 k allocations: 38.457 MiB, 0.17% gc time)
25-th -0.63301548265852 0.007340492709279601
  0.780807 seconds (42.39 k allocations: 38.504 MiB)
26-th -0.6300864856968793 0.006576883266048058
  0.787450 seconds (42.38 k allocations: 38.480 MiB, 0.13% gc time)
27-th -0.6274227436671878 0.005892252532124953
  0.782271 seconds (42.37 k allocations: 38.471 MiB, 0.14% gc time)
28-th -0.6249984477604995 0.005278692772073217
  0.775206 seconds (42.38 k allocations: 38.481 MiB, 0.14% gc time)
29-th -0.6227909814572898 0.004729001813477821
  0.781968 seconds (42.38 k allocations: 38.474 MiB, 0.14% gc time)
30-th -0.6207803405763253 0.004236642241088458
  0.786387 seconds (42.37 k allocations: 38.471 MiB, 0.14% gc time)
31-th -0.6189487129414846 0.0037957051360124226
  0.784794 seconds (42.38 k allocations: 38.479 MiB, 0.12% gc time)
32-th -0.6172801374355674 0.003400858384159586
  0.788380 seconds (42.38 k allocations: 38.498 MiB, 0.12% gc time)
33-th -0.6157602345110169 0.0030473007844825273
  0.796208 seconds (42.38 k allocations: 38.482 MiB)
34-th -0.6143759966380645 0.0027307205386831934
  0.783867 seconds (42.37 k allocations: 38.469 MiB, 0.12% gc time)
35-th -0.6131156094137203 0.0024472404597514195
  0.784699 seconds (42.38 k allocations: 38.487 MiB, 0.12% gc time)
36-th -0.6119683239104643 0.002193392792552506
  0.797728 seconds (42.38 k allocations: 38.476 MiB, 0.12% gc time)
37-th -0.6109243302390619 0.001966063888845468
  0.811322 seconds (42.37 k allocations: 38.466 MiB, 0.12% gc time)
38-th -0.6099746645043274 0.0017624669009970566
  0.798703 seconds (42.38 k allocations: 38.475 MiB, 0.18% gc time)
39-th -0.6091111276751108 0.001580108345611027
  0.804311 seconds (42.38 k allocations: 38.484 MiB, 0.14% gc time)
40-th -0.6083262105712027 0.001416755095863366
  0.802179 seconds (42.37 k allocations: 38.463 MiB, 0.14% gc time)
41-th -0.6076130419061138 0.0012704129923773262
  0.793821 seconds (42.38 k allocations: 38.480 MiB)
42-th -0.6069653199286934 0.00113929572631485
  0.803291 seconds (42.38 k allocations: 38.483 MiB, 0.12% gc time)
43-th -0.6063772737453785 0.0010218056690919832
  0.799485 seconds (42.37 k allocations: 38.476 MiB, 0.14% gc time)
44-th -0.6058436172654165 0.0009165157705025486
  0.797902 seconds (42.38 k allocations: 38.472 MiB, 0.14% gc time)
45-th -0.6053595096728797 0.0008221475370220531
  0.797214 seconds (42.37 k allocations: 38.472 MiB, 0.14% gc time)
46-th -0.6049205193354682 0.0007375584043463334
  0.797498 seconds (42.38 k allocations: 38.481 MiB, 0.13% gc time)
47-th -0.604522590849413 0.0006617273730522861
  0.800408 seconds (42.38 k allocations: 38.492 MiB, 0.12% gc time)
48-th -0.6041620163787921 0.000593740534829054
  0.795930 seconds (42.38 k allocations: 38.491 MiB, 0.11% gc time)
49-th -0.6038354056749607 0.0005327789447378642
  0.794780 seconds (42.37 k allocations: 38.463 MiB, 0.14% gc time)
50-th -0.6035396613195223 0.00047811220217767336
  0.793870 seconds (42.38 k allocations: 38.479 MiB)
51-th -0.603271956206971 0.0004290844742470022
  0.797897 seconds (42.38 k allocations: 38.477 MiB, 0.11% gc time)
52-th -0.6030297098728377 0.0003851109026020435
  0.801658 seconds (42.38 k allocations: 38.478 MiB, 0.13% gc time)
53-th -0.6028105687238502 0.0003456655189107658
  0.797855 seconds (42.38 k allocations: 38.487 MiB, 0.13% gc time)
54-th -0.6026123890946519 0.0003102796361298741
  0.797622 seconds (42.37 k allocations: 38.457 MiB, 0.10% gc time)
55-th -0.6024332167560249 0.00027853279022410063
  0.801957 seconds (42.37 k allocations: 38.468 MiB, 0.13% gc time)
56-th -0.6022712741026642 0.00025004790660750666
  0.801176 seconds (42.37 k allocations: 38.462 MiB, 0.10% gc time)
57-th -0.6021249422723195 0.00022448825972236396
  0.797505 seconds (42.38 k allocations: 38.482 MiB, 0.10% gc time)
58-th -0.6019927502404482 0.00020155103761652473
  0.796861 seconds (42.38 k allocations: 38.479 MiB)
59-th -0.6018733605571956 0.0001809666539966145
  0.803191 seconds (42.38 k allocations: 38.485 MiB, 0.13% gc time)
60-th -0.6017655584431989 0.00016249232065718839
  0.805552 seconds (42.38 k allocations: 38.483 MiB, 0.10% gc time)
61-th -0.6016682408866258 0.00014591021496670423
  0.800285 seconds (42.37 k allocations: 38.467 MiB, 0.11% gc time)
62-th -0.6015804069224595 0.0001310259234340552
  0.804897 seconds (42.38 k allocations: 38.469 MiB, 0.12% gc time)
63-th -0.6015011472832793 0.00011766420738475466
  0.796194 seconds (42.37 k allocations: 38.472 MiB, 0.12% gc time)
64-th -0.6014296388255697 0.00010566954846444353
  0.804151 seconds (42.38 k allocations: 38.481 MiB, 0.12% gc time)
65-th -0.6013651365577349 9.490118655278245e-5
```


