doc"""
solve(solver::Solver, beta::Real, ntau::Integer, ntau_integral::Integer=10ntau)

# Arguments
- `ntau` : number of imaginary time for structure factor S(k,t) and green's function G(r,t)
- `ntau_integral` : number of imaginary time for integration
"""
function solve(solver::Solver, beta::Real, ntau::Integer, ntau_integral::Integer=10*ntau)
    dt = beta/ntau_integral
    L = solver.L
    ef = solver.ef
    nk = length(0:2:L)
    SF = zeros(ntau, nk)
    GF_ca = zeros(ntau, L)
    GF_ac = zeros(ntau, L)
    GK_ca = zeros(ntau, nk)
    GK_ac = zeros(ntau, nk)
    stagN = zeros(nk)
    stagN2 = zeros(nk)
    stagchi = zeros(nk)
    
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

    U1 = similar(rho)
    U2 = similar(rho)
    for (ik, k) in enumerate(0:2:L)
        nstag = ef.vectors' * orderparameter(solver,k) * ef.vectors
        stagN[ik] = tr(nstag*rho)*invZ
        stagN2[ik] = tr(nstag*nstag*rho)*invZ

        for it in 1:ntau_integral
            t1 = dt*(it-1)
            t2 = beta-t1
            U1 .= diagm(0=>exp.(-t1.*ef.values))
            U2 .= diagm(0=>exp.(-t2.*ef.values))
            stagchi[ik] += invZ*tr(U2*nstag*U1*nstag)
        end
    end
    stagchi .*= dt
    stagchi .-= beta.*stagN.^2
    stagchi .*= L

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
            ss = invZ * tr(sf * basis(solver,1))
            GF_ca[it, mod(i-1,L)+1] = -invZ * tr(gfca * annihilator(solver,1))
            GF_ac[it, mod(i-1,L)+1] = -invZ * tr(gfac * creator(solver,1))
            for (ik,k) in enumerate(0:2:L)
                fourier_factor = cospi(k*invV*(i-1))
                SF[it,ik] += fourier_factor * ss
                GK_ca[it,ik] += fourier_factor * GF_ca[it, mod(i-1,L)+1]
                GK_ac[it,ik] += fourier_factor * GF_ac[it, mod(i-1,L)+1]
            end
        end
    end
    SF .-= stagN'.^2

    V2 = L*L
    return Dict("Energy"=>E, "Total Energy"=>E*L,
                "Energy^2"=>E2, "Total Energy^2"=>E2*V2,
                "Specific Heat"=>C, "Heat Capacity"=>C*L,
                "Order Parameter"=>stagN, "Susceptibility"=>stagchi,
                "Structure Factor"=>SF,
                "Temperature Green's Function ca in r"=>GF_ca,
                "Temperature Green's Function ac in r"=>GF_ac,
                "Temperature Green's Function ca in k"=>GK_ca,
                "Temperature Green's Function ac in k"=>GK_ac,
               )
end
