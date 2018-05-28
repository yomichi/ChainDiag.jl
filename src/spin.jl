@inline ldof(S) = Int(2S+1)
function Sp(S::Real)
    S2 = Int(2S)
    ret = [0.5sqrt((S2-m2)*(S2+m2+2) ) for m2 in (S2-2):-2:(-S2)]
    return diagm(ret,1)
end
function Sp(S, L, i)
    leftN = ldof(S)^(i-1)
    rightN = ldof(S)^(L-i)
    kron(eye(leftN), kron(Sp(S), eye(rightN)))
end

function Sm(S::Real)
    S2 = Int(2S)
    ret = [0.5sqrt((S2+m2)*(S2-m2+2) ) for m2 in (S2):-2:(-S2+2)]
    return diagm(ret,-1)
end
function Sm(S, L, i)
    leftN = ldof(S)^(i-1)
    rightN = ldof(S)^(L-i)
    kron(eye(leftN), kron(Sm(S), eye(rightN)))
end

function Sz(S::Real)
    S2 = Int(2S)
    ret = [0.5m2 for m2 in S2:-2:-S2]
    return diagm(ret)
end
function Sz(S, L, i)
    leftN = ldof(S)^(i-1)
    rightN = ldof(S)^(L-i)
    kron(eye(leftN), kron(Sz(S), eye(rightN)))
end

function totalSz(S,L; staggered::Bool=false)
    sz = Sz(S)
    S2 = Int(2S)
    ld = ldof(S)
    N = ld^L
    res = zeros(N)
    for i in 1:N
        for (j,m) in enumerate(digits(i-1,ld,L))
            res[i] += sz[m+1,m+1] * ifelse(staggered && iseven(j), -1.0, 1.0)
        end
    end
    return diagm(res)
end

function spinchain(S, L, Jz, Jxy, h)
    S2 = Int(2S)
    ld = ldof(S)
    sz = Sz(S)
    sp = Sp(S)
    sm = Sm(S)
    jxy = 0.5Jxy
    bondH = Jz.*kron(sz,sz) .+ jxy.*(kron(sp,sm) .+ kron(sm,sp)) - (0.5h).*(kron(sz, eye(ld)) .+ kron(eye(ld), sz))
    N = ld^L
    H = Jz.*kron(sz, kron(eye(div(N,ld*ld)), sz))
    H .+= jxy .* kron(sm, kron(eye(div(N,ld*ld)), sp))
    H .+= jxy .* kron(sp, kron(eye(div(N,ld*ld)), sm))
    H .+= (-0.5h) .* kron(sz, eye(div(N,ld)))
    H .+= (-0.5h) .* kron(eye(div(N,ld)), sz)

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
\mathcal{H} = Jz \sum_i(S^z_i S^z_{i+1}) + 0.5Jxy \sum_i(S^+_i S^-_{i+1} + h.c.) - h \sum_i S^z_i
"""
struct SpinChainSolver
    ef :: Base.LinAlg.Eigen{Float64, Float64, Matrix{Float64}, Vector{Float64}}
    S :: Float64
    L :: Int
    SpinChainSolver(S, L, Jz, Jxy, h) = new(eigfact(spinchain(S, L, Jz, Jxy, h)), S, L)
end


function solve(H::SpinChainSolver, beta::Real, ntau::Integer)
    S = H.S
    L = H.L
    ef = H.ef
    nk = length(0:2:L)
    SS = zeros(nk, ntau)
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
    mz = ef.vectors' * totalSz(S,L) * ef.vectors
    mz2 = mz*mz
    M = trace(mz*rho)*invZ*invV
    M2 = trace(mz2*rho)*invZ*invV^2
    chi = L*beta*(M2-M^2)

    mz .= ef.vectors' * totalSz(S,L, staggered=true) * ef.vectors
    mz2 .= mz*mz
    stagM = trace(mz*rho)*invZ*invV
    stagM2 = trace(mz2*rho)*invZ*invV^2
    stagchi = L*beta*(stagM2-stagM^2)

    for it in 1:ntau
        t1 = beta*((it-1)/ntau)
        t2 = beta-t1
        U1 = ef.vectors * diagm(exp.(-t1.*ef.values)) * ef.vectors'
        U2 = ef.vectors * diagm(exp.(-t2.*ef.values)) * ef.vectors'
        for i in 1:L
            A = U2*Sz(S, L,i)*U1
            for j in 1:L
                ss = invZ * trace(A * Sz(S, L,j))
                for (ik,k) in enumerate(0:2:L)
                    SS[ik,it] += invV * cospi(k*invV*(i-j)) * ss
                end
            end
        end
    end
    V2 = L*L
    return Dict("Energy"=>E, "Total Energy"=>E*L,
                "Energy^2"=>E2, "Total Energy^2"=>E2*V2,
                "Specific Heat"=>C, "Heat Capacity"=>C*L,
                "Magnetization"=>M, "Total Magnetization"=>M*L,
                "Magnetization^2"=>M2, "Total Magnetization^2"=>M2*V2,
                "Susceptibility"=>chi,
                "Staggered Magnetization"=>stagM, "Total Staggered Magnetization"=>stagM*L,
                "Staggered Magnetization^2"=>stagM2, "Total Staggered Magnetization^2"=>stagM2*V2,
                "Staggered Susceptibility"=>stagchi,
                "Structure Factor"=>SS,
               )
end
