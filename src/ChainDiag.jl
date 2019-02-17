__precompile__()

module ChainDiag

using Markdown
using LinearAlgebra

export SpinChainSolver, BosonChainSolver, solve

abstract type Solver end

include("spin.jl")
include("boson.jl")
include("solve.jl")

end # module ChainDiag
