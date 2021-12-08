# Contains structures and functions for simulating decays 
struct DecaySimulationSettings
    # Configuration for simulating decay 
    irf::TcspcTransient
    bck::TcspcTransient
    model::TcspcFitModel
    function DecaySimulationSettings(irf::TcspcTransient,
        bck::TcspcTransient,model::TcspcFitModel)
        # Ensure that the irf and background have the same axis 
        if irf.t != bck.t 
            error("irf and bck must have same time axis")
        end
        # Ensure that the linear model components are nonnegative 
        if !isnonnegative(model)
            error("Linear coefficients must be nonnegative")
        end
        new(irf,bck,model)
    end
end

function sim_component(τ::Real,nphots::Integer,t::TcspcTimeAxis)
    # Simulates single component 

    # Construct Exponential distribution 
    d = Exponential(τ)
    phots = rand(d,nphots)
    # Fit data to histogram 
    bins = time_bins(t)
    h = fit(Histogram,phots,bins)
    return h.weights
end

function sim_noise(nphots::Integer,t::TcspcTimeAxis)
    # Simulates random noise photons 
    bins = time_bins(t)
    d = Uniform(bins[1],bins[end])
    phots = rand(d,nphots)
    h = fit(Histogram,phots,bins)
    return h.weights
end

function simulate_dk(settings::DecaySimulationSettings)
    # Generates a multi-exponential decay w/ IRF convolution

    # First, shift our IRF 
    irf = shift(settings.irf,settings.model.Q.val)
    irf = TcspcTransient(settings.irf.t,irf)
    # Now, determine the number of photons we need per component 
    M = model_vec(settings.model,irf,settings.bck)
    N = sum(M)
    # Multiply N by the fraction of the contribution to determine number of photons 
    # Only the offset and the exponential components need this 
    # To account for scattering and background effects, we have the IRF and background
    # already provided. Hence, we can just multiply the vector by the fraction, and add to the 
    # binned histogram prior to Poisson weighting.
    τₘ = τₛ(settings.model)
    αₑ = τₘ[1:2:end,1]
    τₑ = τₘ[2:2:end,1]
    if nacomp(settings.model) > 0
        θₘ = θₛ(settings.model)
        αₑₘ = θₘ[1:2:end,[1, 4]]
        αₑ = vcat(αₑ,αₑₘ[:,1] .* αₑₘ[:,2])
        τₑₘ = θₘ[2:2:end,[1, 4]]
        τₑ = vcat(τₑ,((1 ./ τₑₘ[:,1]) + (1 ./ τₑₘ[:,2])) .^ (-1))
    end
    ncomp = length(τₑ)
    A = zeros(irf.t.nbins,ncomp)
    for i = 1:ncomp
        τtemp = τₑ[i]
        nphots = round(Int, N * αₑ[i])
        dkᵢ = sim_component(τtemp,nphots,irf.t)
        dkᵢ /= maximum(dkᵢ)
        A[:,i] = dkᵢ
    end
    # Now, convolve matrix in place 
    convmat!(A,irf)
    # Now, multiply by effective weights 
    dk = A[:,1:nlcomp(settings.model)] * αₑ[1:nlcomp(settings.model)]
    if nacomp(settings.model) > 0
        dk += A[:,(nlcomp(settings.model) + 1):end] * 
            αₑ[(nlcomp(settings.model) + 1):end]
    end
    # Now simulate noise 
    nphot = round(Int,settings.model.Z.val * N)
    dk += sim_noise(nphot,irf.t)
    # Now, account for scatter and background 
    dk += irf.data * settings.model.S.val
    dk += settings.bck.data * settings.model.V.val
    # Finally, weight everything using Poisson 
    sim_poisson(round.(Int,dk))
end

function simulate_dk(irf::TcspcTransient,bck::TcspcTransient,
    model::TcspcFitModel)
    # Wrapper for simulation 
    sett = DecaySimulationSettings(irf,bck,model)
    simulate_dk(sett)
end

function simulate_dk(start::AbstractFloat,
    stop::AbstractFloat,nbins::Integer,
    ψ::AbstractFloat, irf::AbstractVector{<:Integer},
    bck::AbstractVector{<:Integer},
    param::AbstractVector{<:AbstractFloat},
    ntau::Integer,nani::Integer)

    # Simulates decays in public API 

    # First, build time axis 
    tax = TcspcTimeAxis(start,stop,nbins,ψ)
    # Next, build irf and bck 
    irf = TcspcTransient(tax,irf)
    bck = TcspcTransient(tax,bck)
    # Finally, build generic model from params
    # Keep all parameters free - no fitting here 
    fxparam = falses(size(param))
    gparam = falses(size(param))
    mod = TcspcFitModel(param,fxparam,gparam,ntau,nani,true)
    simulate_dk(irf,bck,mod)
end

function simulate_dk(tcal::AbstractFloat,nbins::Integer,
    ψ::AbstractFloat, irf::AbstractVector{<:Integer},
    bck::AbstractVector{<:Integer},
    param::AbstractVector{<:AbstractFloat},
    ntau::Integer,nani::Integer)

    # Simulates decays in public API 

    # First, build time axis 
    tax = TcspcTimeAxis(tcal,nbins,ψ)
    # Next, build irf and bck 
    irf = TcspcTransient(tax,irf)
    bck = TcspcTransient(tax,bck)
    # Finally, build generic model from params
    # Keep all parameters free - no fitting here 
    fxparam = falses(size(param))
    gparam = falses(size(param))
    mod = TcspcFitModel(param,fxparam,gparam,ntau,nani,true)
    simulate_dk(irf,bck,mod)
end