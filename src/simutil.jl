# Utility functions for simulating data 
function sim_poisson(trans::AbstractVector{<:Integer})
    # Simulates Poisson noise for a given input vector

    # First, determine the unique distributions from which we need to draw 
    λs = unique(trans)
    N = length(λs)
    for λ ∈ λs
        # Generate different Poisson distribution for each 
        d = Poisson(λ)
        # Find indices to fill 
        ind = findall(==(λ),trans)
        nval = length(ind)
        # Generate nval random number 
        vals = rand(d,nval)
        # Put values back into the array
        trans[ind] = vals
    end
    return trans
end