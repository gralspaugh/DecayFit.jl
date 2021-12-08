# Source file for working with different decays


struct DecayParameter
    # Generic type for configuring sub-module
    val::AbstractFloat 
    fixed::Bool
    glob::Bool 
end

function param_vector(pmain::AbstractVector{<:AbstractFloat},
    fxparam::AbstractVector{<:Bool},
    gparam::AbstractVector{<:Bool})
    
    # Converts guesses to DecayParameters 

    # Ensure sizes match 
    if length(pmain) != length(fxparam) || 
        length(fxparam) != length(gparam)
        error("Parameter vectors must have same length ($(length(pmain))")
    end
    npar = length(pmain)
    param = Vector{DecayParameter}(undef,npar)
    # Fill in vector 
    for i = 1:npar 
        param[i] = DecayParameter(pmain[i],fxparam[i],gparam[i])
    end
    return param
end

struct DecayComponent 
    # Stores amplitude-weighted exponential component 
    α::DecayParameter
    τ::DecayParameter
end

struct TcspcFitModel 
    # Configuration for fitting data 
    τ::AbstractVector{DecayComponent}
    θ::AbstractVector{DecayComponent}
    Q::DecayParameter 
    S::DecayParameter 
    V::DecayParameter
    Z::DecayParameter
    function TcspcFitModel(pmain::AbstractVector{DecayParameter},
        ntau::Integer,nani::Integer)
        # Builds fit model using parameter vector 

        # Ensure that dimensions match specs 
        if length(pmain) != ((2 * ntau) + (2 * nani) + 4)
            error("Parameter vector must be $(((2 * ntau) + (2 * nani) + 4)) elements long")
        end
        if ntau < 1
            error("Must provide at least one lifetime")
        end
        if nani < 0
            error("Number of anisotropies must be nonnegative")
        end
        α = pmain[1:2:(2 * ntau)]
        τ = pmain[2:2:(2 * ntau)]
        τvec = Vector{DecayComponent}(undef,ntau)
        for i = 1:ntau
            τvec[i] = DecayComponent(α[i],τ[i])
        end
        if nani > 0
            sa = (2 * ntau) + 1
            ea = sa + (2 * nani)
            β = pmain[sa:2:ea]
            θ = pmain[(sa + 1):2:ea]
            θvec = Vector{DecayComponent}(undef,nani)
            for i = 1:nani
                θvec[i] = DecayComponent(β[i],θ[i])
            end
        else
            θvec = Vector{DecayComponent}(undef,0)
        end
        Q = pmain[(end - 3)]
        S = pmain[(end - 2)]
        V = pmain[(end - 1)]
        Z = pmain[end]
        new(τvec,θvec,Q,S,V,Z)
    end
end

function TcspcFitModel(pmain::AbstractVector{<:AbstractFloat},
    fxparam::AbstractVector{Bool},
    gparam::AbstractVector{Bool},
    ntau::Integer,nani::Integer)
    # Wrapper for TcspcFitModel construction 
    param = param_vector(pmain,fxparam,gparam)
    TcspcFitModel(param,ntau,nani)
end

function ncomp(mod::TcspcFitModel)
    # Returns number of variable components as tuple 
    τₙ = length(mod.τ)
    θₙ = length(mod.θ)
    return (τₙ,θₙ)
end

function npar(mod::TcspcFitModel)
    # Returns number of parameters in model 
    N = sum(ncomp(mod))
    # Don't include color shift in buffer 
    return (2 * N) + 3
end

function param_buffer(mod::TcspcFitModel)
    # Groups parameters for easy access in fitting 

    pmain = zeros(3,npar(mod))
    # Row 1: Values
    # Row 2: Fixed/free (specified as 1)
    # Row 3: Global/non-global (specified as 1)
    
    # Start with tau 
    τᵢ = 1
    for τ ∈ mod.τ
        pmain[:,2 * τᵢ - 1] = [τ.α.val;τ.α.fixed;τ.α.glob]
        pmain[:,2 * τᵢ] = [τ.τ.val;τ.τ.fixed;τ.τ.glob]
        τᵢ += 1
    end
    if length(mod.θ) > 0
        θᵢ = τᵢ
        for θ ∈ mod.θ
            pmain[:,2 * θᵢ - 1] = [θ.α.val;θ.α.fixed;θ.α.glob]
            pmain[:,2 * θᵢ] = [θ.τ.val;θ.τ.fixed;θ.τ.glob]
            θᵢ += 1
        end
    end
    pmain[:,end - 2] = [mod.S.val;mod.S.fixed;mod.S.glob]
    pmain[:,end - 1] = [mod.V.val;mod.V.fixed;mod.V.glob]
    pmain[:,end] = [mod.Z.val;mod.Z.fixed;mod.Z.glob]
    return transpose(pmain)
end

function isnonnegative(model::TcspcFitModel)
    # Checks whether linear coefficients are nonnegative 
    tf = true
    nlcomp = length(model.τ)
    for τ ∈ model.τ
        if τ.α.val < 0
            tf = false 
            return tf 
        end
    end  
    nacomp = length(model.θ)
    for θ ∈ model.θ
        if θ.α.val < 0
            tf = false 
            return tf 
        end
    end
    tf = tf && model.S.val >= 0 && model.V.val >= 0 && model.Z.val >= 0
end

function taumat(model::TcspcFitModel,tcal::AbstractFloat,
    nbins::Integer)
    # Column for each lifetime in model 
    t = -t_vector(tcal,nbins)
    nlcomp = length(model.τ)
    A = zeros(nbins,nlcomp)
    Ai = 1
    for τ ∈ model.τ
        A[:,Ai] = exp.(t / τ.τ.val)
        Ai += 1
    end
    return A
end

function efftau(model::TcspcFitModel)
    # Returns matrix of lifetimes with amplitudes 

    nlcomp = length(model.τ)
    x = zeros(nlcomp,2)
    xi = 1
    for τ ∈ model.τ
        x[xi,1] = τ.α.val
        x[xi,2] = τ.τ.val
        xi += 1
    end
    return x
end

function effanitau(model::TcspcFitModel)
    # Returns matrix of effective lifetimes with anisotropy components

    nlcomp = length(model.τ)
    nacomp = length(model.θ)
    ncomp = nlcomp * nacomp 
    x = zeros(ncomp,2)
    xi = 1
    for θ ∈ model.θ 
        for τ ∈ model.τ
            x[xi,1] = τ.α.val * θ.α.val
            x[xi,2] = ((1 / τ.τ.val) + (1 / θ.τ.val))
            xi += 1
        end
    end
    x[:,2] = x[:,2] .^ (-1)
    return x
end

function animat(model::TcspcFitModel,tcal::AbstractFloat,
    nbins::Integer)

    t = -t_vector(tcal,nbins)
    τeff = effanitau(model)
    ncomp = size(τeff,1)
    A = zeros(nbins,ncomp)
    for i = 1:ncomp
        A[:,i] = exp.(t / τeff[i,2])
    end
    return A
end

function combmat(model::TcspcFitModel,tcal::AbstractFloat,
    nbins::Integer)
    # Returns the combined tau and anisotropy tau matrix 
    A1 = taumat(model,tcal,nbins)
    if length(model.θ) > 0
        A2 = animat(model,tcal,nbins)
        return [A1, A2]
    else
        return A1
    end
end

function combmat(model::TcspcFitModel,tax::TcspcTimeAxis)
    # Does the exact same thing, except with a time axis wrapper
    combmat(model,tax.tcal,tax.nbins)
end

function convmat(model::TcspcFitModel,irf::TcspcTransient)
    # Returns the convolved anisotropy matrix from model 
    A = combmat(model,irf.t)
    # Convolve each column with the IRF and replace it
    for i = 1:size(A,2)
        temp = conv(irf.data,A[:,i])
        # Try interpolating to see if things change
        A[:,i] = linitpconv(temp,irf.t)
        # A[:,i] = temp[1:irf.t.nbins]
    end
    return A
end

function linmat(model::TcspcFitModel,irf::TcspcTransient,
    bck::TcspcTransient)
    # Returns full unweighted linear matrix 
    if irf.t != bck.t
        error("irf and bck must be same!")
    end
    A = convmat(model,irf)
    return [A, irf.data, bck.data, ones(irf.t.nbins)]
end

function linitpconv(temp::AbstractVector{<:AbstractFloat},
    t::TcspcTimeAxis)
    tnew = range(t.tcal,t.tcal * t.nbins,length=length(temp))
    itp = LinearInterpolation(tnew,temp)
    itp(t_vector(t))
end

function wconvmat(model::TcspcFitModel,irf::TcspcTransient)
    # Returns the weighted convolved anisotropy matrix from model 

    A = convmat(model,irf)
    # First few: pure lifetime 
    Ai = 1
    for τ ∈ model.τ
        A[:,Ai] *= τ.α.val
        Ai += 1
    end
    # Next: combined anisotropy and lifetime 
    for θ ∈ model.θ 
        for τ ∈ model.τ
            A[:,Ai] *= (θ.α.val * τ.α.val)
            Ai += 1
        end
    end
    return A
end

function wconv_vec(model::TcspcFitModel,irf::TcspcTransient)
    # Returns the weighted vector convolution 

    A = wconvmat(model,irf)
    nlcomp = length(model.τ)
    # Start by summing up the pure lifetimes
    I = (1 / 3) * sum(A[:,1:nlcomp],dims=2)
    # Now, do the same except for anisotropies 
    I += (η(irf.t) / 3) * sum(A[:,(nlcomp+1):end],dims=2)
end

function model_vec(model::TcspcFitModel,irf::TcspcTransient,
    bck::TcspcTransient)
    # Returns total model vector 
    if irf.t != bck.t
        error("irf and bck must have same time axis!")
    end
    I = wconv_vec(model,irf)
    A1 = model.S.val * irf.data
    A2 = model.V.val * bck.data
    A3 = model.Z.val
    I .+ A1 .+ A2 .+ A3
end




