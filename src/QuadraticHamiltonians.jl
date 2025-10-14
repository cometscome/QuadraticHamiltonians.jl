module QuadraticHamiltonians
export QuadraticOPs, FermionOP, Hamiltonian,
    construct_matrix

include("fermionoperators.jl")

include("Hamiltonians.jl")

include("./meanfields/RSCGSolver.jl")
export RSCGSolver, solve

include("./meanfields/ChebyshevSolver.jl")
export ChebyshevSolver

include("./meanfields/Meanfields.jl")
export Meanfields_solver

include("./projection/projection.jl")

end
