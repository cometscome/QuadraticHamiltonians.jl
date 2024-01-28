
using RSCG
using SparseIR
import SparseIR: valueim
export RSCGSolver

struct RSCGSolver{T1,T2}
    T::Float64 #temperature
    wmax::Float64
    smpl_beta::T1
    smpl_Matsubara::T2
    ωn_s::Vector{ComplexF64}
end

function RSCGSolver(T; kargs...)#wmax=10.0)
    beta = 1 / T
    if :wmax in keys(kargs)
        wmax = values(kargs).wmax
    else
        wmax = 10.0
    end

    basis = FiniteTempBasis(Fermionic(), beta, wmax, 1e-7)
    smpl = MatsubaraSampling(basis)
    ωn_s = valueim.(smpl.sampling_points, beta)
    println("num. of Matsubara freqs. ", length(ωn_s))
    smpl_beta = TauSampling(basis; sampling_points=[beta])
    return RSCGSolver{typeof(smpl_beta),typeof(smpl)}(T, wmax, smpl_beta, smpl, ωn_s)
end

function fit_ir(rscg::RSCGSolver, Gij)
    return fit_ir(Gij, rscg.smpl_Matsubara, rscg.smpl_beta)
end

function fit_ir(Gij, smpl_Matsubara, smpl_beta)
    gl = fit(smpl_Matsubara, Gij)
    G0 = evaluate(smpl_beta, gl)
    return -G0[1]
end

function solve(rscg::RSCGSolver, A, i, j)
    #println(rscg.ωn_s)
    Gij = greensfunctions(i, j, rscg.ωn_s, A)
    Gij0 = fit_ir(Gij, rscg.smpl_Matsubara, rscg.smpl_beta)
    return Gij0
end
