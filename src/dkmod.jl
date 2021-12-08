# Contains code helpful for using transient decay models 

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
    ntau::Integer,nani::Integer,Q::Bool)
    # Wrapper for TcspcFitModel construction 
    param = param_vector(pmain,fxparam,gparam)
    if !Q
        # Add a temporary, fixed, non-global color shift with val = 0
        q = DecayParameter(0,true,false)
        np = length(param) + 1
        ptemp = Vector{DecayParameter}(undef,np)
        ind = vcat(1:(np - 4),(np - 2):np)
        ptemp[ind] = param
        ptemp[(np - 3)] = q
        param = ptemp
    end
    TcspcFitModel(param,ntau,nani)
end

nlcomp(mod::TcspcFitModel) = length(mod.τ)

nacomp(mod::TcspcFitModel) = length(mod.θ)

npar(mod::TcspcFitModel) = 2 * (nlcomp(mod) + nacomp(mod)) + 3

function npar(arr::AbstractVector{TcspcFitModel})
    # Vectorized npar 
    np = zeros(Int,length(arr))
    modcount = 1
    for mod ∈ arr
        np[modcount] = npar(mod)
        modcount += 1
    end
    return np
end

function param_buffer(mod::TcspcFitModel)
    # Groups parameters for easy access in fitting 

    pmain = zeros(npar(mod),3)
    # Column 1: Values
    # Column 2: Fixed/free (specified as 1)
    # Column 3: Global/non-global (specified as 1)
    
    # Start with tau 
    τᵢ = 1
    for τ ∈ mod.τ
        pmain[2 * τᵢ - 1,:] = [τ.α.val,τ.α.fixed,τ.α.glob]
        pmain[2 * τᵢ,:] = [τ.τ.val,τ.τ.fixed,τ.τ.glob]
        τᵢ += 1
    end
    if nacomp(mod) > 0
        θᵢ = τᵢ
        for θ ∈ mod.θ
            pmain[2 * θᵢ - 1,:] = [θ.α.val,θ.α.fixed,θ.α.glob]
            pmain[2 * θᵢ,:] = [θ.τ.val,θ.τ.fixed,θ.τ.glob]
            θᵢ += 1
        end
    end
    pmain[end - 2,:] = [mod.S.val,mod.S.fixed,mod.S.glob]
    pmain[end - 1,:] = [mod.V.val,mod.V.fixed,mod.V.glob]
    pmain[end,:] = [mod.Z.val,mod.Z.fixed,mod.Z.glob]
    return pmain
end

function isunique(arr::AbstractVector{TcspcFitModel})
    # Determines whether the models are the same 

    # Let's assume they are at first 
    for i = 2:length(arr)
        modₐ = arr[i - 1]
        modᵦ = arr[i]
        # First, check if length of parameters is the same 
        if ncomp(modₐ) != ncomp(modᵦ)
            return false
        end
        # Now, check lifetimes for consistency 
        for j = 1:length(modₐ.τ)
            # Amplitude check 
            if modₐ.τ[j].α.fixed != modᵦ.τ[j].α.fixed ||
                modₐ.τ[j].α.glob != modᵦ.τ[j].α.glob
                return false 
            end
            # Lifetime check 
            if modₐ.τ[j].τ.fixed != modᵦ.τ[j].τ.fixed ||
                modₐ.τ[j].τ.glob != modᵦ.τ[j].τ.glob
                return false 
            end
        end
        # Now, check anisotropies for consistency 
        for j = 1:length(modₐ.θ)
            # Amplitude check 
            if modₐ.θ[j].α.fixed != modᵦ.θ[j].α.fixed ||
                modₐ.θ[j].α.glob != modᵦ.θ[j].α.glob
                return false 
            end
            # Lifetime check 
            if modₐ.θ[j].τ.fixed != modᵦ.θ[j].τ.fixed ||
                modₐ.θ[j].τ.glob != modᵦ.θ[j].τ.glob
                return false 
            end
        end
        # Finally, check the others 
        if modₐ.Q.glob != modᵦ.Q.glob
            return false
        end
        if modₐ.S.fixed != modᵦ.S.fixed ||
            modₐ.S.glob != modᵦ.S.glob
            return false 
        end
        if modₐ.V.fixed != modᵦ.V.fixed ||
            modₐ.V.glob != modᵦ.V.glob
            return false 
        end
        if modₐ.Z.fixed != modᵦ.Z.fixed ||
            modₐ.Z.glob != modᵦ.Z.glob
            return false 
        end
    end
    return true
end

function isnonnegative(model::TcspcFitModel)
    # Checks whether linear coefficients are nonnegative 
    tf = true
    for τ ∈ model.τ
        if τ.α.val < 0
            tf = false 
            return tf 
        end
    end  
    for θ ∈ model.θ
        if θ.α.val < 0
            tf = false 
            return tf 
        end
    end
    tf = tf && model.S.val >= 0 && model.V.val >= 0 && model.Z.val >= 0
end
