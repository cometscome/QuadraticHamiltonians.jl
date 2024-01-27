# QuadraticHamiltonians

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://cometscome.github.io/QuadraticHamiltonians.jl/stable/)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://cometscome.github.io/QuadraticHamiltonians.jl/dev/)
[![Build Status](https://github.com/cometscome/QuadraticHamiltonians.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/cometscome/QuadraticHamiltonians.jl/actions/workflows/CI.yml?query=branch%3Amain)


This is an alpha version. 

```julia
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
H = -1.0C_{1,1}C_{2,1}^+ -1.0C_{1,1}C_{16,1}^+ -1.0C_{2,1}C_{3,1}^+ -1.0C_{2,1}C_{1,1}^+ -1.0C_{3,1}C_{4,1}^+ -1.0C_{3,1}C_{2,1}^+ -1.0C_{4,1}C_{5,1}^+ -1.0C_{4,1}C_{3,1}^+ -1.0C_{5,1}C_{6,1}^+ -1.0C_{5,1}C_{4,1}^+ -1.0C_{6,1}C_{7,1}^+ -1.0C_{6,1}C_{5,1}^+ -1.0C_{7,1}C_{8,1}^+ -1.0C_{7,1}C_{6,1}^+ -1.0C_{8,1}C_{9,1}^+ -1.0C_{8,1}C_{7,1}^+ -1.0C_{9,1}C_{10,1}^+ -1.0C_{9,1}C_{8,1}^+ -1.0C_{10,1}C_{11,1}^+ -1.0C_{10,1}C_{9,1}^+ -1.0C_{11,1}C_{12,1}^+ -1.0C_{11,1}C_{10,1}^+ -1.0C_{12,1}C_{13,1}^+ -1.0C_{12,1}C_{11,1}^+ -1.0C_{13,1}C_{14,1}^+ -1.0C_{13,1}C_{12,1}^+ -1.0C_{14,1}C_{15,1}^+ -1.0C_{14,1}C_{13,1}^+ -1.0C_{15,1}C_{16,1}^+ -1.0C_{15,1}C_{14,1}^+ -1.0C_{16,1}C_{1,1}^+ -1.0C_{16,1}C_{15,1}^+ 
---------------------------------
16×16 SparseArrays.SparseMatrixCSC{Float64, Int64} with 32 stored entries:
⎡⠪⡢⡀⠀⠀⠀⠀⠈⎤
⎢⠀⠈⠪⡢⡀⠀⠀⠀⎥
⎢⠀⠀⠀⠈⠪⡢⡀⠀⎥
⎣⡀⠀⠀⠀⠀⠈⠪⡢⎦
```