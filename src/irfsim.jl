# Contains functions for simulating IRFs and exponential decays

struct IrfSimulationSettings
    # Configuration for simulating IRF
    t::TcspcTimeAxis
    N::Unsigned
    fwhm::AbstractFloat
    offset::AbstractFloat
    # Inner constructor for argument validation
    function IrfSimulationSettings(t::TcspcTimeAxis,
        N::Integer,fwhm::AbstractFloat,
        offset::AbstractFloat)
        if fwhm < 0
            error("fwhm is less than 0")
        end
        if offset < 0
            error("offset is less than 0")
        end
        new(t,N,fwhm,offset)
    end
end

function simulate_irf(settings::IrfSimulationSettings)
    # Generates a Gaussian instrument response function 

    σ = settings.fwhm/(2 * sqrt(2 * log(2)))
    gaussDist = Normal(settings.offset, σ)
    # Generate photons
    photons = rand(gaussDist,settings.N)
    # Fit histogram to photons at time points
    bins = time_bins(settings.t)
    h = fit(Histogram,photons,bins)
    # Weights contain the number of counts per bin, which is what we want
    return sim_poisson(h.weights)
end

function simulate_irf(tax::TcspcTimeAxis,
    N::Integer,fwhm::AbstractFloat,
    offs::AbstractFloat)
    sett = IrfSimulationSettings(tax,N,fwhm,offs)
    simulate_irf(sett)
end

function simulate_irf(start::AbstractFloat,
    stop::AbstractFloat,nbins::Integer,
    ψ::AbstractFloat,
    N::Integer,fwhm::AbstractFloat,
    offs::AbstractFloat)
    tax = TcspcTimeAxis(start,stop,nbins,ψ)
    simulate_irf(tax,N,fwhm,offs)
end

function simulate_irf(tcal::AbstractFloat,nbins::Integer,
    ψ::AbstractFloat, N::Integer,fwhm::AbstractFloat,
    offs::AbstractFloat)
    # Open method for dispatch 
    tax = TcspcTimeAxis(tcal,nbins,ψ)
    simulate_irf(tax,N,fwhm,offs)
end