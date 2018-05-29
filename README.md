# ChainDiag.jl

Fulldiag solver for finite temperature problem of XXZ spin model and softcore Bose Hubbard model on a periodic chain lattice.

Spin length `S` and the maximum number of partcles on each site `M` can take any value, provided sufficiently large RAM and fast CPU.

## Install
```julia
julia> Pkg.clone("https://github.com/yomichi/ChainDiag.jl")
```

## Short sample

The following example solves hardcore Bose-Hubbard model:

``` julia
julia> using ChainDiag

julia> M=1; L=6; t=1.0; V=0.0; U=0.0; mu=0.0;

# hardcore Bose-Hubbard model
# H = -\sum_{i=1}^6 [c_i a_{i+1} + a_i c_{i+1}]
julia> solver = BosonChainSolver(M, L, t=t, V=V, U=U, mu=mu);
# default parameters: t=1.0, V=U=mu=0.0

julia> beta=10.0; ntau=10;

julia> res = solve(solver, 10.0, 10); # => Dict{String,Any}

julia> res["Energy"] # energy density
-0.6658339097820527

julia> res["Total Energy"]
-3.995003458692316

julia> res["Specific Heat"]
0.04421188790391106

julia> res["Order Parameter"] # particle density
0.4999999999999998

julia> res["Susceptibility"] # for order parameter
0.0155393999976694

julia> res["Staggered Order Parameter"]
-1.142555719539883e-16

julia> res["Staggered Susceptibility"]
4.984460600002333

# structure factor for order parameter
# S[ik,itau], where k = (ik-1)*pi/L and tau = beta*(itau-1)/ntau
julia> res["Structure Factor"]
4×10 Array{Float64,2}:
 1.50155   1.50155    1.50155      1.50155      1.50155      1.50155     1.50155      1.50155      1.50155      1.50155
 0.166667  0.0226205  0.00307281   0.000418025  5.79739e-5   1.55311e-5  5.79739e-5   0.000418025  0.00307281   0.0226205
 0.333333  0.0167645  0.000868712  4.94162e-5   3.607e-6     7.42208e-7  3.607e-6     4.94162e-5   0.000868712  0.0167645
 0.498446  0.0478126  0.00610676   0.000819932  0.000112829  2.99855e-5  0.000112829  0.000819932  0.00610676   0.0478126

julia> res["Temperature Green's Function ac in r"] # G[ir,itau] = - <a(ir,tau) c(1,0)>
6×10 Array{Float64,2}:
 -0.5       -0.194659  -0.107583  -0.0666032  -0.0467005  -0.040697   -0.0467005  -0.0666032  -0.107583  -0.194659
 -0.332917  -0.180486  -0.106152  -0.0664507  -0.0466832  -0.0406929  -0.0466832  -0.0664507  -0.106152  -0.180486
 -0.222222  -0.167477  -0.104765  -0.0663039  -0.0466672  -0.0406894  -0.0466672  -0.0663039  -0.104765  -0.167477
 -0.221945  -0.167342  -0.104727  -0.0662943  -0.0466648  -0.0406883  -0.0466648  -0.0662943  -0.104727  -0.167342
 -0.222222  -0.167477  -0.104765  -0.0663039  -0.0466672  -0.0406894  -0.0466672  -0.0663039  -0.104765  -0.167477
 -0.332917  -0.180486  -0.106152  -0.0664507  -0.0466832  -0.0406929  -0.0466832  -0.0664507  -0.106152  -0.180486

julia> res["Temperature Green's Function ac in k"] # G[ik,itau]
4×10 Array{Float64,2}:
 -1.83222    -1.05793     -0.634144    -0.398407     -0.280066    -0.24415     -0.280066    -0.398407     -0.634144    -1.05793
 -0.38875    -0.0403257   -0.0042435   -0.000455742  -5.17996e-5  -1.23097e-5  -5.17996e-5  -0.000455742  -0.0042435   -0.0403257
 -0.166805   -0.0140386   -0.00139326  -0.000142987  -1.49222e-5  -3.05381e-6  -1.49222e-5  -0.000142987  -0.00139326  -0.0140386
 -0.0566659  -0.00129847  -8.31761e-5  -1.52699e-5   -3.63198e-6  -1.59243e-6  -3.63198e-6  -1.52699e-5   -8.31761e-5  -0.00129847
```

To solve spin chain, use `SpinChainSolver(S, L; Jz, Jxy, h)` instead of `BoseChainSolver(M, L; t, V, U, mu)`.

