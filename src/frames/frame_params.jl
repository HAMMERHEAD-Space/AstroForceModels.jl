# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#
# Description
# ==========================================================================================
#
#   Parameter types for frame-aware force models
#
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

export FrameAwareParams

"""
    FrameAwareParams{F,E}

Parameter struct for frame-aware force model propagation.

This struct wraps around a standard `ComponentVector` and adds frame transformation
capabilities. It's designed to work seamlessly with both new frame-aware code and
legacy code that doesn't use frames.

# Fields
- `params::ComponentVector`: Standard parameter vector (contains JD and other numeric params)
- `frames::F`: FrameSystem for coordinate transformations
- `epoch::E`: Reference epoch (Tempo.Epoch)
- `propagation_frame::Symbol`: Frame in which the state vector is expressed

# Usage
```julia
using Tempo, FrameTransformations

# Create frames
frames = create_default_frames()

# Create epoch
epoch = Epoch("2023-01-01T00:00:00 TDB")

# Create base params
base_params = ComponentVector(JD = 2460000.5)

# Wrap with frame info
params = FrameAwareParams(base_params, frames, epoch, :ICRF)

# Access like normal ComponentVector
params.JD  # Works!

# Access frame info
params.frames
params.epoch
params.propagation_frame
```

# Notes
- Can be indexed like a ComponentVector (delegates to the inner params)
- Property access first checks frame fields, then delegates to params
- Fully compatible with ODE solvers
"""
struct FrameAwareParams{F,E,P<:ComponentVector}
    params::P
    frames::F
    epoch::E
    propagation_frame::Symbol
end

# Delegation to inner ComponentVector
Base.getindex(p::FrameAwareParams, i) = getindex(p.params, i)
Base.setindex!(p::FrameAwareParams, v, i) = setindex!(p.params, v, i)
Base.length(p::FrameAwareParams) = length(p.params)
Base.size(p::FrameAwareParams) = size(p.params)
Base.iterate(p::FrameAwareParams, args...) = iterate(p.params, args...)

# Property access: check frame fields first, then delegate to params
function Base.getproperty(p::FrameAwareParams, s::Symbol)
    if s in (:params, :frames, :epoch, :propagation_frame)
        return getfield(p, s)
    else
        return getproperty(getfield(p, :params), s)
    end
end

function Base.setproperty!(p::FrameAwareParams, s::Symbol, v)
    if s in (:params, :frames, :epoch, :propagation_frame)
        error("Cannot modify frame-aware parameter fields directly")
    else
        return setproperty!(getfield(p, :params), s, v)
    end
end

Base.propertynames(p::FrameAwareParams) = tuple(:frames, :epoch, :propagation_frame, propertynames(p.params)...)
Base.eltype(p::FrameAwareParams) = eltype(p.params)
Base.keys(p::FrameAwareParams) = keys(p.params)
Base.values(p::FrameAwareParams) = values(p.params)
Base.pairs(p::FrameAwareParams) = pairs(p.params)
Base.similar(p::FrameAwareParams) = FrameAwareParams(similar(p.params), p.frames, p.epoch, p.propagation_frame)
Base.similar(p::FrameAwareParams, ::Type{T}) where {T} = FrameAwareParams(similar(p.params, T), p.frames, p.epoch, p.propagation_frame)
Base.axes(p::FrameAwareParams) = axes(p.params)
Base.IndexStyle(::Type{<:FrameAwareParams}) = IndexLinear()

# Show method
function Base.show(io::IO, p::FrameAwareParams)
    println(io, "FrameAwareParams:")
    println(io, "  propagation_frame: ", p.propagation_frame)
    println(io, "  epoch: ", p.epoch)
    println(io, "  frames: ", typeof(p.frames))
    println(io, "  params: ", p.params)
end

