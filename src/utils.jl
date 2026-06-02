"""
    current_jd(p, t)

Compute the current Julian Date from parameters `p` and elapsed time `t` [s].
"""
@inline current_jd(p, t) = p.JD + t / 86400.0

"""
    angle_between_vectors(
        v1::AbstractVector{T1}, v2::AbstractVector{T2}
    ) where {T1<:Number,T2<:Number}

Computes the angle between two vectors in a more numerically stable way than dot product.

# Arguments
-`v1::AbstractVector{<:Number}`: The first vector of the computation
-`v2::AbstractVector{<:Number}`: The second vector of the computation

# Returns 
-`angle::Number`: The angle between the two vectors
"""
@inline function angle_between_vectors(
    v1::AbstractVector{T1}, v2::AbstractVector{T2}
) where {T1<:Number,T2<:Number}
    T = promote_type(T1, T2)

    unitv1 = v1 ./ _euclidean_norm(v1)
    unitv2 = v2 ./ _euclidean_norm(v2)

    y = unitv1 - unitv2
    x = unitv1 + unitv2

    a = 2.0 * atan(_euclidean_norm(y), _euclidean_norm(x))

    angle::T = if !(sign(a) == -1.0 || sign(T(π) - a) == -1.0)
        a
    else
        (sign(a) == -1.0 ? zero(T) : T(π))
    end

    return angle
end

# Euclidean (2-) norm computed with an explicit reduction. This intentionally
# avoids `LinearAlgebra.norm`/`normalize`, whose `StridedVector{<:Union{...}}`
# methods route through `ReinterpretArray` padding code that hits a Base bug on
# Julia 1.12 when analyzed against the generic `AbstractVector{<:Number}`
# signature used here.
@inline function _euclidean_norm(v::AbstractVector{T}) where {T<:Number}
    s = abs2(zero(T))
    @inbounds for i in eachindex(v)
        s += abs2(v[i])
    end
    return sqrt(s)
end
