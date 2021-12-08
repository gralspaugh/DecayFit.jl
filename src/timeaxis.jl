# Contains types for time axis manipulations 
@def_structequal struct TcspcTimeAxis
    # Stores time bin information for the collected axes
    nbins::Unsigned
    tcal::AbstractFloat
    ψ::AbstractFloat
    # Define internal constructor 
    function TcspcTimeAxis(tcal::AbstractFloat,nbins::Integer,
        ψ::AbstractFloat)
        # Round nbins to nearest power of two
        if nbins <= 0
            error("nbins ($nbins) must be greater than 0")
        end
        nbins = nextpow(2,nbins)
        if tcal < 0
            error("tcal must be larger than 0.0")
        end
        new(nbins,tcal,ψ)
    end
end

# Multiple dispatch with inputs of start, stop, and nbins 
function TcspcTimeAxis(start::AbstractFloat,
    stop::AbstractFloat,nbins::Integer,
    ψ::AbstractFloat)
    if start >= stop
        error("start ($start) must be less than stop ($stop)")
    end
    # Determine tcal 
    tcal = (stop - start) / nbins
    TcspcTimeAxis(tcal,nbins,ψ)
end

# Helper methods with time axis 
function chan2time(tax::TcspcTimeAxis,chan::Integer)
    # Converts channel number to time point along axis
    if chan <= 0
        error("chan must be greater than 0")
    end
    # Return to middle of bin
    chan * (tax.tcal / 2)
end

function time2chan(tax::TcspcTimeAxis,t::AbstractFloat)
    if t < 0
        error("t must be greater than or equal to 0")
    end
    # Generate histogram of data
    h = fit(Histogram,t,time_bins(tax))
    # Find non-zero elements
    findfirst(h.weights)
end

time_bins(t::TcspcTimeAxis) = range(0,step=t.tcal,length=t.nbins+1)

t_vector(t::TcspcTimeAxis) = range(t.tcal,step=t.tcal,length=t.nbins)

function t_vector(tcal::AbstractFloat,nbins::Integer)
    tax = TcspcTimeAxis(tcal,nbins,44.7) # placeholder ψ
    range(tax.tcal,step=tax.tcal,length=nbins)
end

η(t::TcspcTimeAxis) = 3 * (cosd(t.ψ) ^ 2) - 1

