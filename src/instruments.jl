# Functions for working with irfs 

function shift(data::AbstractVector{<:Real},Q::Integer)
    # Shifts the irf by n channels
    nbins = length(data)
    temp = data
    if Q > 0 
        nsc = Q
        nec = nbins
        osc = 1
        oec = nbins - Q
        temp[1:(nsc - 1)] .= 0
    elseif Q < 0
        nsc = 1
        nec = nbins - Q 
        osc = abs(Q)
        oec = nbins
        temp[nec:nbins] .= 0
    end
    temp[nsc:nec] = data[osc:oec]
end

function shift(irf::TcspcTransient,Q::AbstractFloat)
    if Q != 0
        # Same thing, except rounds Q to nearest tenth 
        Q = round(Int,round(Q,digits=1) / 0.1)
        # Linear IRF interpolation
        t = t_vector(irf.t)
        itp = interpolate((t,),irf.data,Gridded(Linear()))
        # Now evaluate at the new grid
        tnew = TcspcTimeAxis(irf.t.tcal / 10,irf.t.nbins * 10,irf.t.ψ)
        tgrid = t_vector(tnew)
        v = itp(tgrid)
        v = shift(v,Q)
        # Now, reverse interpolate back 
        itp = interpolate((tgrid,),v,Gridded(Linear()))
        # Evaluate at the old grid
        round(Int,itp(t))
    else
        irf.data
    end
    
end


