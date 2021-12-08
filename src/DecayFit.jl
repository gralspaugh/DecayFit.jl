module DecayFit

# Add external packages here

using Random, Distributions, StatsBase
using StructEquality, FFTW

# Add file dependencies here
include("timeaxis.jl")
include("transients.jl")
include("instruments.jl")
include("dkmod.jl")
include("dkmat.jl")
include("simutil.jl")
include("irfsim.jl")
include("dksim.jl")

# Add function exports here
export simulate_irf, simulate_dk

# Write your package code here. 

end
