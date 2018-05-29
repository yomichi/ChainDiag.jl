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
    GF_ca = zeros(L, ntau)
    GF_ac = zeros(L, ntau)
    GK_ca = zeros(nk, ntau)
    GK_ac = zeros(nk, ntau)
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
    N = trace(n*rho)*invZ
    N2 = trace(n2*rho)*invZ
    chi = L*beta*(N2-N^2)

    n = ef.vectors' * orderparameter(solver,true) * ef.vectors
    n2 = n*n
    stagN = trace(n*rho)*invZ
    stagN2 = trace(n2*rho)*invZ
    stagchi = L*beta*(stagN2-stagN^2)

    for it in 1:ntau
        t1 = beta*((it-1)/ntau)
        t2 = beta-t1
        U1 = ef.vectors * diagm(exp.(-t1.*ef.values)) * ef.vectors'
        U2 = ef.vectors * diagm(exp.(-t2.*ef.values)) * ef.vectors'
        for i in 1:L
            sf = U2*basis(solver,i)*U1
            gfca = U2*creator(solver,i)*U1
            gfac = U2*annihilator(solver,i)*U1
            ss = invZ * trace(sf * basis(solver,1))
            GF_ca[mod(i-1,L)+1, it] = -invZ * trace(gfca * annihilator(solver,1))
            GF_ac[mod(i-1,L)+1, it] = -invZ * trace(gfac * creator(solver,1))
            for (ik,k) in enumerate(0:2:L)
                fourier_factor = cospi(k*invV*(i-1))
                SF[ik,it] += fourier_factor * ss
                GK_ca[ik,it] += fourier_factor * GF_ca[mod(i-1,L)+1, it]
                GK_ac[ik,it] += fourier_factor * GF_ac[mod(i-1,L)+1, it]
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
                "Structure Factor"=>SF,
                "Temperature Green's Function ca in r"=>GF_ca,
                "Temperature Green's Function ac in r"=>GF_ac,
                "Temperature Green's Function ca in k"=>GK_ca,
                "Temperature Green's Function ac in k"=>GK_ac,
               )
end

end # module ChainDiag
