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
    FrameAwareParams{F,E,P}

Parameter struct for frame-aware force model propagation.

Wraps a `FrameSystem`, epoch, propagation frame, and an optional `ComponentVector`
for user-defined parameters. The Julian Date is derived from `epoch` automatically.

# Fields
- `params::ComponentVector`: User-defined parameters (can be empty).
- `frames::F`: FrameSystem for coordinate transformations.
- `epoch::E`: Reference epoch (Tempo.Epoch). JD is derived from this.
- `propagation_frame::Symbol`: Frame in which the state vector is expressed.

# Constructors
```julia
# Minimal — just frames, epoch, and propagation frame
p = FrameAwareParams(frames, epoch, :ICRF)

# With extra user parameters
p = FrameAwareParams(ComponentVector(; custom=1.0), frames, epoch, :ICRF)
```

# Property access
- `p.frames`, `p.epoch`, `p.propagation_frame` → struct fields
- `p.JD` → computed from epoch (read-only, never stale)
- `p.custom` → delegates to inner ComponentVector
"""
struct FrameAwareParams{F,E,P<:ComponentVector}
    params::P
    frames::F
    epoch::E
    propagation_frame::Symbol
end

# Convenience: no ComponentVector needed
function FrameAwareParams(frames::F, epoch::E, propagation_frame::Symbol) where {F,E}
    return FrameAwareParams(ComponentVector(), frames, epoch, propagation_frame)
end

# Delegation to inner ComponentVector
Base.getindex(p::FrameAwareParams, i) = getindex(p.params, i)
Base.setindex!(p::FrameAwareParams, v, i) = setindex!(p.params, v, i)
Base.length(p::FrameAwareParams) = length(p.params)
Base.size(p::FrameAwareParams) = size(p.params)
Base.iterate(p::FrameAwareParams, args...) = iterate(p.params, args...)

# Property access: struct fields and derived JD first, then delegate to params
function Base.getproperty(p::FrameAwareParams, s::Symbol)
    if s in (:params, :frames, :epoch, :propagation_frame)
        return getfield(p, s)
    elseif s === :JD
        # Derive JD from epoch — single source of truth, never stale
        return 2451545.0 + j2000s(getfield(p, :epoch)) / 86400.0
    else
        return getproperty(getfield(p, :params), s)
    end
end

function Base.setproperty!(p::FrameAwareParams, s::Symbol, v)
    if s in (:params, :frames, :epoch, :propagation_frame, :JD)
        error("Cannot modify FrameAwareParams fields directly")
    else
        return setproperty!(getfield(p, :params), s, v)
    end
end

function Base.propertynames(p::FrameAwareParams)
    tuple(:frames, :epoch, :propagation_frame, :JD, propertynames(p.params)...)
end
Base.eltype(p::FrameAwareParams) = Float64
Base.keys(p::FrameAwareParams) = keys(p.params)
Base.values(p::FrameAwareParams) = values(p.params)
Base.pairs(p::FrameAwareParams) = pairs(p.params)
function Base.similar(p::FrameAwareParams)
    FrameAwareParams(similar(p.params), p.frames, p.epoch, p.propagation_frame)
end
function Base.similar(p::FrameAwareParams, ::Type{T}) where {T}
    FrameAwareParams(similar(p.params, T), p.frames, p.epoch, p.propagation_frame)
end
Base.axes(p::FrameAwareParams) = axes(p.params)
Base.IndexStyle(::Type{<:FrameAwareParams}) = IndexLinear()

# Show method
function Base.show(io::IO, p::FrameAwareParams)
    println(io, "FrameAwareParams:")
    println(io, "  propagation_frame: ", p.propagation_frame)
    println(io, "  epoch: ", getfield(p, :epoch))
    println(io, "  JD: ", p.JD, " (derived from epoch)")
    println(io, "  frames: ", typeof(getfield(p, :frames)))
    isempty(p.params) || println(io, "  params: ", p.params)
end
