
import .Mod_RSCGSolver: RSCGSolver
struct Meanfields_solver{M}
    method::M
end

function Meanfields_solver(ham::Hamiltonian, T, ; method="RSCG", kargs...)
    if method == "RSCG"
        rscg = RSCGSolver(T; kargs...)
        return Meanfields_solver{typeof(rscg)}(rscg)
    else
        error("method $method is not supported yet")
    end
end
