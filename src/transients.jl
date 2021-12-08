# Contains structures and methods related to transient types

struct TcspcTransient
    data::AbstractVector{UInt32}
    t::TcspcTimeAxis
    # Inner constructor for main validation
    function TcspcTransient(tax::TcspcTimeAxis,
        data::AbstractVector{<:Integer})
        # Ensure that number of bins match 
        nbins = length(data)
        if nbins != tax.nbins
            error("Number of elements ($nbins) must match length of axis ($(tax.nbins))")
        end
        new(data,tax)
    end
end

function TcspcTransient(tcal::AbstractFloat,
    data::AbstractVector{<:Integer},
    ψ::AbstractFloat)
    # Determine number of bins from data 
    nbins = length(data)
    if nbins != nextpow(2,nbins)
        error("nbins ($nbins) must be power of 2")
    end
    # Build time axis 
    t = TcspcTimeAxis(tcal,nbins,ψ)
    TcspcTransient(t,data)
end

# Multiple dispatch for padded input 
function TcspcTransient(start::AbstractFloat,
    stop::AbstractFloat,data::AbstractVector{<:Integer},
    ψ::AbstractFloat)
    # Construct time axis from start, stop, 
    # number of data elements
    nbins = length(data)
    tax = TcspcTimeAxis(start,stop,nbins,ψ)
    # Determine whether we need to pad zeros 
    if tax.nbins != nbins
        # Determine where to start padding on both sides
        leftPadEnd = time2chan(tax,start)
        rightPadStart = time2chan(tax,stop)
        # Construct a new data vector
        temp = zeros(eltype(data),size(tax))
        # Fill new array with data
        temp[leftPadEnd:(rightPadStart - 1)] = data
        data = temp 
    end
    TcspcTransient(tax.tcal,data)
end

# Helper methods 
function poisson_weights(trans::TcspcTransient)
    # Returns 1 / data for all counts above 15
    data = zeros(size(trans.data))
    data[findall(<=(15),trans.data)] .= 1/15
    data[findall(==(0),trans.data)] .= 1
    dind = findall(>(15),trans.data)
    data[dind] = 1 ./ trans.data[dind]
    return data
end

function guess_start_stop(trans::TcspcTransient)
    # Guesses where the start of the curve is, 
    # as well as where to stop fitting data
    maxtup = findmax(trans.data)
    # Find point of highest rise in curve
    x = 1:trans.t.nbins
    y = trans.data
    itp = interpolate((x,),y,Gridded(Linear()))
    deriv = only.(Interpolations.gradient.(Ref(itp),x))
    dmaxtup = findmax(deriv)
    # Good approx: 2 * dmaxtup[2] - maxtup[2]
    fit_start = 2 * dmaxtup[2] - maxtup[2]
    # Find first element less than 50 
    fit_end = findfirst(<(50),y)
    (fit_start,fit_end)
end

function estimate_irf(trans::TcspcTransient)
    # Estimates instrument response function for curve
    maxtup = findmax(trans.data)
    chantup = guess_start_stop(trans)
    # Offset ≈ tcal * fit_start + (maxtup[2] - fit_start) / 2
    offschan = round(Int,chantup[1] + (maxtup[2] - chantup[1]) / 2)
    # Determine time
    offs = chan2time(trans.t,offschan)
    # FWHM ≈ tcal * (maxtup[2] - fit_start)
    fwhm = chan2time(trans.t,round(Int,maxtup[2] - chantup[1]))
    # N = sum(data)
    N = sum(trans.data)
    # Simulate IRF, then rescale it
    irf = simulate_irf(trans.t,N,fwhm,offs)
    mxirftup = findmax(irf)
    # Scale = mxIRF / mxData
    sc = mxirftup[1] / maxtup[1]
    round(Int,sc * irf)
end

