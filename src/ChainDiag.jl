__precompile__()

module ChainDiag

export SpinChainSolver, BosonChainSolver, solve

abstract type Solver end

include("spin.jl")
include("boson.jl")
include("solve.jl")

end # module ChainDiag
