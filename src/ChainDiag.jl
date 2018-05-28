__precompile__()

module ChainDiag

export SpinChainSolver, BosonChainSolver, solve

abstract type Solver end

include("spin.jl")
include("boson.jl")

function solve(solver::Solver, beta::Real, ntau::Integer)
    L = solver.L
    ef = solver.ef
    nk = length(0:2:L)
    SF = zeros(nk, ntau)
    CF = zeros(L, ntau)
    Z = 0.0
    E = 0.0
    E2 = 0.0
    invV = 1.0/L
    for en in reverse(ef.values)
        z = exp(-beta*en)
        Z += z
        E += z*(en*invV)
        E2 += z*(en*invV)^2
    end
    invZ = 1.0/Z
    E *= invZ
    E2 *= invZ
    C = L*beta^2*(E2-E^2)

    rho = diagm(exp.(-beta.*ef.values))
    n = ef.vectors' * orderparameter(solver) * ef.vectors
    n2 = n*n
    N = trace(n*rho)*invZ*invV
    N2 = trace(n2*rho)*invZ*invV^2
    chi = L*beta*(N2-N^2)

    n = ef.vectors' * orderparameter(solver,true) * ef.vectors
    n2 = n*n
    stagN = trace(n*rho)*invZ*invV
    stagN2 = trace(n2*rho)*invZ*invV^2
    stagchi = L*beta*(stagN2-stagN^2)

    for it in 1:ntau
        t1 = beta*((it-1)/ntau)
        t2 = beta-t1
        U1 = ef.vectors * diagm(exp.(-t1.*ef.values)) * ef.vectors'
        U2 = ef.vectors * diagm(exp.(-t2.*ef.values)) * ef.vectors'
        for i in 1:L
            A = U2*basis(solver,i)*U1
            B = U2*creator(solver,i)*U1
            for j in 1:L
                ss = invZ * trace(A * basis(solver,j))
                for (ik,k) in enumerate(0:2:L)
                    SF[ik,it] += invV * cospi(k*invV*(i-j)) * ss
                end
                CF[mod(i-j,L)+1, it] += invZ * trace(B * annihilator(solver,j))
            end
        end
    end
    V2 = L*L
    return Dict("Energy"=>E, "Total Energy"=>E*L,
                "Energy^2"=>E2, "Total Energy^2"=>E2*V2,
                "Specific Heat"=>C, "Heat Capacity"=>C*L,
                "Order Parameter"=>N, "Order Parameter^2"=>N2,
                "Susceptibility"=>chi,
                "Staggered Order Parameter"=>stagN, "Staggered Order Parameter^2"=>stagN2,
                "Staggered Susceptibility"=>stagchi,
                "Structure Factor"=>SF, "Correlation Function"=>CF,
               )
end

end # module ChainDiag
