doc"""
solve(solver::Solver, beta::Real, ntau::Integer, ntau_integral::Integer=10ntau)

# Arguments
- `ntau` : number of imaginary time for structure factor S(k,t) and green's function G(r,t)
- `ntau_integral` : number of imaginary time for integration
"""
function solve(solver::Solver, beta::Real, ntau::Integer, ntau_integral::Integer=10*ntau)
    L = solver.L
    ef = solver.ef
    nk = length(0:2:L)
    SF = zeros(ntau, nk)
    GF_ca = zeros(ntau, L)
    GF_ac = zeros(ntau, L)
    GK_ca = zeros(ntau, nk)
    GK_ac = zeros(ntau, nk)
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

    rho = diagm(0=>exp.(-beta.*ef.values))
    n = ef.vectors' * orderparameter(solver) * ef.vectors
    n2 = n*n
    N = trace(n*rho)*invZ
    N2 = trace(n2*rho)*invZ
    chi = L*beta*(N2-N^2)

    nstag = ef.vectors' * orderparameter(solver,true) * ef.vectors
    stagN = trace(nstag*rho)*invZ
    stagN2 = trace(nstag*nstag*rho)*invZ

    chi = 0.0
    stagchi = 0.0
    dt = beta/ntau_integral
    U1 = similar(n)
    U2 = similar(n)
    for it in 1:ntau_integral
        t1 = dt*(it-1)
        t2 = beta-t1
        U1 .= diagm(0=>exp.(-t1.*ef.values))
        U2 .= diagm(0=>exp.(-t2.*ef.values))
        chi += invZ*trace(U2*n*U1*n)
        stagchi += invZ*trace(U2*nstag*U1*nstag)
    end
    chi *= dt
    chi -= beta*N^2
    chi *= L
    stagchi *= dt
    stagchi -= beta*stagN^2
    stagchi *= L

    ## structure factor and green's function
    for it in 1:ntau
        t1 = beta*((it-1)/ntau)
        t2 = beta-t1
        U1 .= ef.vectors * diagm(0=>exp.(-t1.*ef.values)) * ef.vectors'
        U2 .= ef.vectors * diagm(0=>exp.(-t2.*ef.values)) * ef.vectors'
        for i in 1:L
            sf = U2*basis(solver,i)*U1
            gfca = U2*creator(solver,i)*U1
            gfac = U2*annihilator(solver,i)*U1
            ss = invZ * trace(sf * basis(solver,1))
            GF_ca[it, mod(i-1,L)+1] = -invZ * trace(gfca * annihilator(solver,1))
            GF_ac[it, mod(i-1,L)+1] = -invZ * trace(gfac * creator(solver,1))
            for (ik,k) in enumerate(0:2:L)
                fourier_factor = cospi(k*invV*(i-1))
                SF[it,ik] += fourier_factor * ss
                GK_ca[it,ik] += fourier_factor * GF_ca[it, mod(i-1,L)+1]
                GK_ac[it,ik] += fourier_factor * GF_ac[it, mod(i-1,L)+1]
            end
        end
    end

    V2 = L*L
    return Dict("Energy"=>E, "Total Energy"=>E*L,
                "Energy^2"=>E2, "Total Energy^2"=>E2*V2,
                "Specific Heat"=>C, "Heat Capacity"=>C*L,
                "Order Parameter"=>N, "Susceptibility"=>chi,
                "Staggered Order Parameter"=>stagN, "Staggered Susceptibility"=>stagchi,
                "Structure Factor"=>SF,
                "Temperature Green's Function ca in r"=>GF_ca,
                "Temperature Green's Function ac in r"=>GF_ac,
                "Temperature Green's Function ca in k"=>GK_ca,
                "Temperature Green's Function ac in k"=>GK_ac,
               )
end
