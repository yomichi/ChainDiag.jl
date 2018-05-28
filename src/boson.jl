# @inline ldof(S) = Int(2S+1)
function creator(M::Integer)
    ret = Float64[sqrt(m) for m in M:-1:1]
    return diagm(ret,1)
end
function creator(M::Integer, L::Integer, i::Integer)
    ld = M+1
    leftN = ld^(i-1)
    rightN = ld^(L-i)
    kron(eye(leftN), kron(creator(M), eye(rightN)))
end

function annihilator(M::Integer)
    ret = Float64[sqrt(m+1) for m in M-1:-1:0]
    return diagm(ret,-1)
end
function annihilator(M::Integer, L::Integer, i::Integer)
    ld = M+1
    leftN = ld^(i-1)
    rightN = ld^(L-i)
    kron(eye(leftN), kron(annihilator(M), eye(rightN)))
end

function number(M::Integer)
    ret = Float64[m for m in M:-1:0]
    return diagm(ret)
end
function number(M::Integer, L::Integer, i::Integer)
    ld = M+1
    leftN = ld^(i-1)
    rightN = ld^(L-i)
    kron(eye(leftN), kron(number(M), eye(rightN)))
end

function totalnumber(M::Integer,L::Integer; staggered::Bool=false)
    ld = M+1
    N = ld^L
    n = number(M)
    res = zeros(N)
    for i in 1:N
        for (j,m) in enumerate(digits(i-1,ld,L))
            res[i] += n[m+1,m+1] * ifelse(staggered && iseven(j), -1.0, 1.0)
        end
    end
    return diagm(res)
end

function bosonchain(M::Integer, L::Integer, t::Real, V::Real, U::Real, mu::Real)
    ld = M+1
    c = creator(M)
    a = annihilator(M)
    n = number(M)
    i = eye(ld)

    onsite = (-0.5mu).*n .+ (0.5U).*n.*(n.-1.0)
    bondH = V.*kron(n,n) .- t.*(kron(c,a) .+ kron(a,c)) .+ kron(onsite,i) .+ kron(i,onsite)

    N = ld^L
    H = V.*kron(n, kron(eye(div(N,ld*ld)), n))
    H .-= t .* kron(c, kron(eye(div(N,ld*ld)), a))
    H .-= t .* kron(a, kron(eye(div(N,ld*ld)), c))
    H .+= kron(onsite, eye(div(N,ld)))
    H .+= kron(eye(div(N,ld)), onsite)

    leftN = 1
    rightN = ld^(L-2)
    for i in 1:L-1
        H .+= kron(eye(leftN), kron(bondH, eye(rightN)))
        leftN *= ld
        rightN = div(rightN,ld)
    end
    return H
end

doc"""
\mathcal{H} = -t\sum_i(a_i c_{i+1} + c_i a_{i+1}) + V\sum_i n_i n_{i+1} + U \sum_i n_i(n_i-1) - mu \sum_i n_i
"""
struct BosonChainSolver
    ef :: Base.LinAlg.Eigen{Float64, Float64, Matrix{Float64}, Vector{Float64}}
    M :: Int
    L :: Int
    BosonChainSolver(M::Integer, L::Integer, t::Real, V::Real, U::Real, mu::Real) = new(eigfact(bosonchain(M, L, t, V, U, mu)), M, L)
end

function solve(H::BosonChainSolver, beta::Real, ntau::Integer)
    M = H.M
    L = H.L
    ef = H.ef
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
    n = ef.vectors' * totalnumber(M,L) * ef.vectors
    n2 = n*n
    N = trace(n*rho)*invZ*invV
    N2 = trace(n2*rho)*invZ*invV^2
    chi = L*beta*(N2-N^2)

    n = ef.vectors' * totalnumber(M,L, staggered=true) * ef.vectors
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
            A = U2*number(M, L,i)*U1
            B = U2*creator(M, L,i)*U1
            for j in 1:L
                ss = invZ * trace(A * number(M, L,j))
                for (ik,k) in enumerate(0:2:L)
                    SF[ik,it] += invV * cospi(k*invV*(i-j)) * ss
                end
                CF[mod(i-j,L)+1, it] += invZ * trace(B * annihilator(M, L,j))
            end
        end
    end
    V2 = L*L
    return Dict("Energy"=>E, "Total Energy"=>E*L,
                "Energy^2"=>E2, "Total Energy^2"=>E2*V2,
                "Specific Heat"=>C, "Heat Capacity"=>C*L,
                "Number Density"=>N, "Number of Particles"=>N*L,
                "Number Density^2"=>N2, "Number of Particles^2"=>N2*V2,
                "Susceptibility"=>chi,
                "Staggered Number Density"=>stagN, "Staggered Number of Particles"=>stagN*L,
                "Staggered Number Density^2"=>stagN2, "Staggered Number of Particles^2"=>stagN2*V2,
                "Staggered Susceptibility"=>stagchi,
                "Structure Factor"=>SF, "Correlation Function"=>CF,
               )
end
