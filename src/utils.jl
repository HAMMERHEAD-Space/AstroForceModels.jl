"""
    current_jd(p, t)

Compute the current Julian Date from parameters `p` and elapsed time `t` [s].
"""
@inline current_jd(p, t) = p.JD + t / 86400.0

@inline function angle_between_vectors(
    v1::AbstractVector{T1}, v2::AbstractVector{T2}
) where {T1,T2}
    T = promote_type(T1, T2)

    unitv1 = SVector{length(v1),T1}(normalize(v1))
    unitv2 = SVector{length(v2),T2}(normalize(v2))

    y = unitv1 - unitv2
    x = unitv1 + unitv2

    a = 2.0 * atan(norm(y), norm(x))

    angle::T =
        !(signbit(a) || signbit(float(T)(π) - a)) ? a : (signbit(a) ? zero(T) : float(T)(π))

    return angle
end
