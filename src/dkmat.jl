# Contains code helpful for the generation of transients 

function τₛ(model::TcspcFitModel)
    # Returns lifetime values with amplitudes as a (2 * N) x 3 array 

    τₘ = zeros(2 * nlcomp(model),3)
    τₛ!(τₘ,model)
    return τₘ
end

function τₛ!(τₘ::AbstractMatrix{<:AbstractFloat},
    model::TcspcFitModel)
    τᵢ = 1
    for τ ∈ model.τ
        τₘ[τᵢ,:] = [τ.α.val,τ.α.fixed,τ.α.glob]
        τₘ[τᵢ + 1,:] = [τ.τ.val,τ.τ.fixed,τ.τ.glob]
        τᵢ += 2
    end
end

function θₛ(model::TcspcFitModel)
    # Returns anisotropy values with amplitudes as a (2 * N * M) x 6 array

    θₘ = zeros(2 * (nacomp(model) * nlcomp(model)),6)
    θₛ!(θₘ,model)
    return θₘ
end

function θₛ!(θₘ::AbstractMatrix{<:AbstractFloat},
    model::TcspcFitModel)
    θᵢ = 1
    for θ ∈ model.θ
        for τ ∈ model.τ
            θₘ[θᵢ,:] = vcat([τ.α.val,τ.α.fixed,τ.α.glob],
                [θ.α.val,θ.α.fixed,θ.α.glob])
            θₘ[θᵢ + 1,:] = vcat([τ.τ.val,τ.τ.fixed,τ.τ.glob],
                [θ.τ.val,θ.τ.fixed,θ.τ.glob])
            θᵢ += 2
        end
    end
end

function taumat(model::TcspcFitModel,tcal::AbstractFloat,
    nbins::Integer)
    # Column for each lifetime in model 
    t = -t_vector(tcal,nbins)
    nlcomp = length(model.τ)
    A = zeros(nbins,nlcomp)
    τₘ = τₛ(model)
    taumat!(A,t,vec(τₘ[2:2:end,1]))
    return A
end

function w_taumat(model::TcspcFitModel,tcal::AbstractFloat,
    nbins::Integer)
    # Weights the matrix by the amplitude 
    t = -t_vector(tcal,nbins)
    nlcomp = length(model.τ)
    A = zeros(nbins,nlcomp)
    τₘ = τₛ(model)
    α = τₘ[1:2:end,1]
    τ = τₘ[2:2:end,1]
    w_taumat!(A,t,α,τ)
    return A
end

function taumat!(A::AbstractVecOrMat{<:AbstractFloat},
    t::AbstractVector{<:AbstractFloat},
    τ::AbstractVector{<:AbstractFloat})
    # Evaluates in place for everything 
    Ai = 1
    for τᵢ ∈ τ
        A[:,Ai] = exp.(t / τᵢ)
        Ai += 1
    end
end

function taumat!(model::TcspcFitModel,
    τₘ::AbstractMatrix{<:AbstractFloat},
    A::AbstractVecOrMat{<:AbstractFloat},
    t::AbstractVector{<:AbstractFloat})
    # Evaluates in place with buffer inputs 
    τₛ!(τₘ,model)
    taumat!(A,t,τₘ[2:2:end,1])
end

function w_taumat!(A::AbstractVecOrMat{<:AbstractFloat},
    t::AbstractVector{<:AbstractFloat},
    α::AbstractVector{<:AbstractFloat},
    τ::AbstractVector{<:AbstractFloat})
    # Evaluates weighted tau matrix in place 
    taumat!(A,t,τ)
    A *= α'
end

function w_taumat!(model::TcspcFitModel,
    τₘ::AbstractMatrix{<:AbstractFloat},
    A::AbstractVecOrMat{<:AbstractFloat},
    t::AbstractVector{<:AbstractFloat})
    τₛ!(τₘ,model)
    taumat!(A,t,τₘ[1:2:end,1],τₘ[2:2:end,1])
end

function animat(model::TcspcFitModel,tcal::AbstractFloat,
    nbins::Integer)

    t = -t_vector(tcal,nbins)
    # Determine effective lifetimes from anisotropies 
    θₘ = θₛ(model)
    τₑₘ = θₘ[2:2:end,[1, 4]]
    τₑ = ((1 ./ τₑₘ[:,1]) + (1 ./ τₑₘ[:,2])) .^ (-1)
    ncomp = length(τₑ)
    A = zeros(nbins,ncomp)
    taumat!(A,t,τₑ)
    return A
end

function animat!(model::TcspcFitModel,
    θₘ::AbstractMatrix{<:AbstractFloat},
    A::AbstractVecOrMat{<:AbstractFloat},
    t::AbstractVector{<:AbstractFloat})
    θₛ!(θₘ,model)
    τₑₘ = θₘ[2:2:end,[1, 4]]
    τₑ = ((1 ./ τₑₘ[:,1]) + (1 ./ τₑₘ[:,2])) .^ (-1)
    taumat!(A,t,τₑ)
end

function w_animat(model::TcspcFitModel,tcal::AbstractFloat,
    nbins::Integer)
    # Weights the matrix by the amplitude 
    t = -t_vector(tcal,nbins)
    θₘ = θₛ(model)
    αₑₘ = θₘ[1:2:end,[1, 4]]
    αₑ = αₑₘ[:,1] .* αₑₘ[:,2]
    τₑₘ = θₘ[2:2:end,[1, 4]]
    τₑ = ((1 ./ τₑₘ[:,1]) + (1 ./ τₑₘ[:,2])) .^ (-1)
    ncomp = length(τₑ)
    A = zeros(nbins,ncomp)
    w_taumat!(A,t,αₑ,τₑ)
    return A
end

function w_animat!(model::TcspcFitModel,
    θₘ::AbstractMatrix{<:AbstractFloat},
    A::AbstractVecOrMat{<:AbstractFloat},
    t::AbstractVector{<:AbstractFloat})
    θₛ!(θₘ,model)
    αₑₘ = θₘ[1:2:end,[1, 4]]
    αₑ = αₑₘ[:,1] .* αₑₘ[:,2]
    τₑₘ = θₘ[2:2:end,[1, 4]]
    τₑ = ((1 ./ τₑₘ[:,1]) + (1 ./ τₑₘ[:,2])) .^ (-1)
    w_taumat!(A,t,αₑ,τₑ)
end

function convmat!(A::AbstractVecOrMat,irf::TcspcTransient)
    # Convolves matrix with instrument response function 
    for i = 1:size(A,2)
        A[:,i] = abs.((ifft(fftshift(fft(A[:,i]) .* fft(irf.data)))))
    end
end

function model_vec(model::TcspcFitModel,irf::TcspcTransient,
    bck::TcspcTransient)
    # Returns model vector summed up
    A = taumat(model,irf.t.tcal,irf.t.nbins)
    # Convolve with irf
    convmat!(A,irf)
    τₘ = τₛ(model)
    m = A * τₘ[1:2:end,1]
    if nacomp(model) > 0
        m *= (1 / 3)
        Aₐ = animat(model,irf.t.tcal,irf.t.nbins)
        convmat!(Aₐ,irf)
        θₘ = θₛ(model)
        αₑₘ = θₘ[1:2:end,[1, 4]]
        αₑ = αₑₘ[:,1] .* αₑₘ[:,2]
        m += (η(irf.t) / 3) * (Aₐ * αₑ)
    end
    m += model.S.val * irf.data
    m += model.V.val * bck.data
    m .+= model.Z.val
    return m
end

function model_vec_gs(model::TcspcFitModel,irf::TcspcTransient,
    bck::TcspcTransient)
    # Returns model vector summed up
    A = taumat(model,irf.t.tcal,irf.t.nbins)
    fill!(A,0)
    # Convolve with grinvald-steinberg
    τₘ = τₛ(model)
    grin_stein_conv!(A,irf.t.tcal,irf,τₘ[1:2:end,1],τₘ[2:2:end,1])
    m = sum(A,dims=2)
    if nacomp(model) > 0
        m *= (1 / 3)
        Aₐ = animat(model,irf.t.tcal,irf.t.nbins)
        fill!(Aₐ,0)
        θₘ = θₛ(model)
        αₑₘ = θₘ[1:2:end,[1, 4]]
        αₑ = αₑₘ[:,1] .* αₑₘ[:,2]
        τₑₘ = θₘ[2:2:end,[1, 4]]
        τₑ = ((1 ./ τₑₘ[:,1]) + (1 ./ τₑₘ[:,2])) .^ (-1)
        grin_stein_conv!(Aₐ,irf.t.tcal,irf,αₑ,τₑ)
        m += (η(irf.t) / 3) * sum(Aₐ,dims=2)
    end
    m += model.S.val * irf.data
    m += model.V.val * bck.data
    m .+= model.Z.val
    return m
end

function model_vec!(model::TcspcFitModel,
    m::AbstractVector{<:AbstractFloat},
    τₘ::AbstractMatrix{<:AbstractFloat},
    A::AbstractVecOrMat{<:AbstractFloat},
    t::AbstractVector{<:AbstractFloat},
    irf::TcspcTransient,bck::TcspcTransient)
    # Evaluates everything in place with anisotropy
    τₛ!(τₘ,model)
    taumat!(A,t,τₘ[2:2:end,1])
    convmat!(A,irf)
    m[:] = A * τₘ[1:2:end,1]
    m += model.S.val * irf.data
    m += model.V.val * bck.data
    m += model.Z.val
end

function model_vec!(model::TcspcFitModel,
    m::AbstractVector{<:AbstractFloat},
    τₘ::AbstractMatrix{<:AbstractFloat},
    θₘ::AbstractMatrix{<:AbstractFloat},
    A::AbstractVecOrMat{<:AbstractFloat},
    t::AbstractVector{<:AbstractFloat},
    irf::TcspcTransient,bck::TcspcTransient)
    # Evaluates everything in place 
    τₛ!(τₘ,model)
    taumat!(A[:,1:nlcomp(model)],t,τₘ[2:2:end,1])
    convmat!(A[:,1:nlcomp(model)],irf)
    m[:] = (1 / 3) * A[:,1:nlcomp(model)] * τₘ[1:2:end,1]
    θₛ!(θₘ,model)
    animat!(model,θₘ,A[:,(nlcomp(model) + 1):end],t)
    convmat!(A[:,(nlcomp(model) + 1):end],irf)
    αₑₘ = θₘ[1:2:end,[1, 4]]
    αₑ = αₑₘ[:,1] .* αₑₘ[:,2]
    m += (η(irf.t) / 3) * (A[:,(nlcomp(model) + 1):end] * αₑ)
    m += model.S.val * irf.data
    m += model.V.val * bck.data
    m += model.Z.val
end