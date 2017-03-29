module transform

export rot2, rotx, roty, rotz
export trot2, trotx, troty, trotz
export se2, se3
export r2t, t2r
export rpy2r, rpy2t, rpy2jac, tr2rpy

typealias RealArray{T <: Real} Array{T,1}
typealias RealMatrix{T <: Real} Matrix{T}
typealias RealMatrix3D{T <: Real} Array{T,3}


"""
    check_argument_units(units)

Assert `units` is valid.
"""
function check_argument_units(units::Symbol)
    @assert (units in (:rad, :deg)) "Expected rad or deg for units, got $units"
end

"""
    check_argument_axis_order(axis_order)

Assert `axis_order` is valid.
"""
function check_argument_axis_order(axis_order::Symbol)
    @assert (axis_order in (:xyz, :zyx)) "Expected xyz or zyx for axis_order, got $axis_order"
end

"""
    check_argument_axis(axis)

Assert `axis` is valid.
"""
function check_argument_axis(axis::Symbol)
    @assert (axis in (:x, :y, :z)) "Expected x, y, or z for axis, got $axis"
end

"""
    convert_angle(θ, units)

Convert `θ` to degrees if `units` is `:deg`, otherwise return `θ`.
"""
function convert_angle(θ::Union{Real, RealArray}, units::Symbol)
    check_argument_units(units)

    if units == :deg
            θ = deg2rad(θ)
    end

    return θ
end


# Rotation matrix generation
# -----------------------------------------------------------------------------

"""
    rot2(θ, [units])

Generate 2 x 2 rotation matrix representing a rotation of `θ`.

Interpret `θ` in radians unless `units` is `:deg`.
"""
function rot2(θ::Real, units::Symbol=:rad)
    θ = convert_angle(θ, units)

    cos_θ = cos(θ)
    sin_θ = sin(θ)

    return [cos_θ -sin_θ;
            sin_θ cos_θ]
end

"""
    rot_any(θ, axis, [units])

Generate 3 x 3 rotation matrix representing a rotation of `θ` about `axis`.

Interpret `θ` in radians unless `units` is `:deg`.
"""
function rot_any(θ::Real, axis::Symbol, units::Symbol=:rad)
    check_argument_axis(axis)

    θ = convert_angle(θ, units)

    cos_θ = cos(θ)
    sin_θ = sin(θ)

    if axis == :x
        return [1 0 0;
                0 cos_θ -sin_θ;
                0 sin_θ cos_θ]
    elseif axis == :y
        return [cos_θ 0 sin_θ;
                0 1 0;
                -sin_θ 0 cos_θ]
    else
        return [cos_θ -sin_θ 0;
                sin_θ cos_θ 0;
                0 0 1]
    end
end

"""
    rotx(θ, [units])

Generate 3 x 3 rotation matrix representing a rotation of `θ` about the x-axis.

Interpret `θ` in radians unless `units` is `:deg`.
"""
rotx(θ::Real, units::Symbol=:rad) = rot_any(θ, :x, units)

"""
    roty(θ, [units])

Generate 3 x 3 rotation matrix representing a rotation of `θ` about the y-axis.

Interpret `θ` in radians unless `units` is `:deg`.
"""
roty(θ::Real, units::Symbol=:rad) = rot_any(θ, :y, units)

"""
    rotz(θ, [units])

Generate 3 x 3 rotation matrix representing a rotation of `θ` about the z-axis.

Interpret `θ` in radians unless `units` is `:deg`.
"""
rotz(θ::Real, units::Symbol=:rad) = rot_any(θ, :z, units)


# Homogeneous transform generation
# -----------------------------------------------------------------------------

"""
    r2t(rot_mat)

Generate a homogenous transform from rotation matrix `rot_mat`.

`rot_mat` can be either 2 x 2 or 3 x 3.
- If `rot_mat` is 2 x 2, the result is 3 x 3.
- If `rot_mat` is 3 x 3, the result is 4 x 4.

If `rot_mat` is a rotation matrix sequence then a homogeneous transform sequence
is returned.

The resulting transform has a zero translational component.
"""
function r2t(rot_mat::Union{RealMatrix,RealMatrix3D})
    if ndims(rot_mat) == 3
        return cat(3, [r2t(rot_mat[:, :, i]) for i in 1:size(rot_mat, 3)]...)
    end

    @assert (size(rot_mat) in ((2, 2), (3, 3))) "Expected array of size (2, 2) or (3, 3), instead had size $(size(rot_mat))."

    is_2d = size(rot_mat) == (2, 2)
    if is_2d
        return [rot_mat [0; 0]; [0 0 1]]
    else
        return [rot_mat [0; 0; 0]; [0 0 0 1]]
    end
end

"""
    t2r(trans_mat)

Generate a rotation matrix from homogeneous transform `trans_mat`.

`trans_mat` can be either 3 x 3 or 4 x 4.
- If `trans_mat` is 3 x 3, the result is 2 x 2.
- If `trans_mat` is 4 x 4, the result is 3 x 3.

If `trans_mat` is a homogeneous transform sequence then a rotation matrix sequence
is returned.

The validity of the rotational part is not checked.
"""
function t2r(trans_mat::Union{RealMatrix, RealMatrix3D})
    if ndims(trans_mat) == 3
        return cat(3, [t2r(trans_mat[:, :, i]) for i in 1:size(trans_mat, 3)]...)
    end

    @assert (size(trans_mat) in ((3, 3), (4, 4))) "Expected array of size (3, 3) or (4, 4), instead had size $(size(trans_mat))."

    return trans_mat[1:end - 1, 1:end - 1]
end

"""
    trot2(θ, [units])

Generate a 3 x 3 homogeneous transform representing a rotation of `θ`.

The translational component is zero.

Interpret `θ` in radians unless `units` is `:deg`.
"""
trot2(θ::Real, units::Symbol=:rad) = r2t(rot2(θ, units))

"""
    trot_any(θ, axis, [units])

Generate 4 x 4 homogeneous transform representing a rotation of `θ` about `axis`.

The translational component is zero.

Interpret `θ` in radians unless `units` is `:deg`.
"""
function trot_any(θ::Real, axis::Symbol, units::Symbol=:rad)
    check_argument_axis(axis)

    rot_func = if axis == :x
                   rotx
               elseif axis == :y
                   roty
               else
                   rotz
               end

    return r2t(rot_func(θ, units))
end

"""
    trotx(θ, [units])

Generate a 4 x 4 homogeneous transform representing a rotation of `θ` about the x-axis.

The translational component is zero.

Interpret `θ` in radians unless `units` is `:deg`.
"""
trotx(θ::Real, units::Symbol=:rad) = trot_any(θ, :x, units)

"""
    troty(θ, [units])

Generate a 4 x 4 homogeneous transform representing a rotation of `θ` about the y-axis.

The translational component is zero.

Interpret `θ` in radians unless `units` is `:deg`.
"""
troty(θ::Real, units::Symbol=:rad) = trot_any(θ, :y, units)

"""
    trotz(θ, [units])

Generate a 4 x 4 homogeneous transform representing a rotation of `θ` about the z-axis.

The translational component is zero.

Interpret `θ` in radians unless `units` is `:deg`.
"""
trotz(θ::Real, units::Symbol=:rad) = trot_any(θ, :z, units)

"""
    se2(x, y, [θ, [units]])

Generate planar transformation and rotation information.

Return a 3 x 3 homogeneous transform representing translation `x` and `y` and
rotation `θ` in the plane.

Interpret `θ` in radians unless `units` is `:deg`.
"""
function se2(x::Real, y::Real, θ::Real=0, units::Symbol=:rad)
    θ = convert_angle(θ, units)

    sin_θ = sin(θ)
    cos_θ = cos(θ)

    return [cos_θ -sin_θ x;
            sin_θ cos_θ y;
            0 0 1]
end

"""
    se3(trans_mat)

Lift 3 x 3 transform to 4 x 4.

Return a 4 x 4 homogeneous transform that represents the same x, y translation and
z rotation as does `trans_mat` (3 x 3).
"""
function se3(trans_mat::RealMatrix)
    @assert (size(trans_mat) == (3, 3)) "Expected array of size (3, 3), instead had size $(size(trans_mat))."

    return [trans_mat[1:2, 1:2] [0; 0] trans_mat[1:2, 3];
            0 0 1 0;
            0 0 0 1]
end


# Conversion between roll/pitch/yaw and rotation matrices/homogeneous transforms
# ------------------------------------------------------------------------------

"""
    rpy2r(roll, pitch, yaw, [units, [axis_order]])

Roll-pitch-yaw angles to rotation matrix.

Return a 3 x 3 orthonormal rotation matrix equivalent to the specified `roll`,
`pitch`, and `yaw` angles.  These correspond to rotations about the x, y, and
z axes respectively.

Interpret `roll`, `pitch`, and `yaw` in radians unless `units` is `:deg`.

Return solutions for sequential rotations about x, y, z axes unless `axis_order`
is `:zyx`, in which case return solutions for sequential rotations about z, y,
x axes.
"""
function rpy2r(roll::Real, pitch::Real, yaw::Real, units::Symbol=:rad, axis_order::Symbol=:xyz)
    check_argument_axis_order(axis_order)

    if axis_order == :xyz
        rot_func_a, rot_func_b, rot_func_c = (rotx, roty, rotz)
    elseif axis_order == :zyx
        rot_func_a, rot_func_b, rot_func_c = (rotz, roty, rotx)
    end

    return rot_func_a(roll, units) * rot_func_b(pitch, units) * rot_func_c(yaw, units)
end

"""
    rpy2t(roll, pitch, yaw, [units, [axis_order]])

Roll-pitch-yaw angles to homogeneous transform.

Return a 4 x 4 homogeneous transform equivalent to the specified `roll`, `pitch`,
and `yaw` angles.  These correspond to rotations about the x, y, and z axes
respectively.

Interpret `roll`, `pitch`, and `yaw` in radians unless `units` is `:deg`.

Return solutions for sequential rotations about x, y, z axes unless `axis_order`
is `:zyx`, in which case return solutions for sequential rotations about z, y,
x axes.
"""
rpy2t(roll::Real, pitch::Real, yaw::Real, units::Symbol=:rad, axis_order::Symbol=:xyz) = r2t(rpy2r(roll, pitch, yaw, units, axis_order))

"""
    rpy2jac(roll, pitch, yaw, [units])

Generate 3 x 3 Jacobian matrix from roll-pitch-yaw angle rates to angular velocity.

Used in the creation of an analytical Jacobian.

Interpret `roll`, `pitch`, and `yaw` in radians unless `units` is `:deg`.
"""
function rpy2jac(roll::Real, pitch::Real, yaw::Real, units::Symbol=:rad)
    roll, pitch, yaw = convert_angle([roll, pitch, yaw], units)

    return [1 0 sin(pitch);
            0 cos(roll) (-cos(pitch) * sin(roll));
            0 sin(roll) (cos(pitch) * cos(roll))]
end

"""
    tr2rpy(mat, [units, [axis_order]])

Generate roll-pitch-yaw angles from homogeneous transform or rotation matrix.

If `mat` is a sequence then a sequence of roll-pitch-yaw angles is returned.

There is a singularity for the case where pitch = π / 2 in which case roll is
aribtrarily set to zero and yaw is the sum (roll + yaw).

Note that textbooks (Paul, Spong) use the rotation order zyx.

Return `roll`, `pitch`, and `yaw` in radians unless `units` is `:deg`.

Return solutions for sequential rotations about x, y, z axes unless `axis_order`
is `:zyx`, in which case return solutions for sequential rotations about z, y,
x axes.
"""
function tr2rpy(mat::Union{RealMatrix, RealMatrix3D}, units::Symbol=:rad, axis_order::Symbol=:xyz)
    if ndims(mat) == 3
        return cat(3, [tr2rpy(mat[:, :, i]) for i in 1:size(mat, 3)]...)
    end

    @assert (size(mat) in ((3, 3), (4, 4))) "Expected array of size (3, 3) or (4, 4), instead had size $(size(mat))."

    if axis_order == :xyz
        singularity_present = all([abs(elem) < eps() for elem in [mat[2, 3] mat[3, 3]]])
        if singularity_present
            roll = 0
            pitch = atan2(mat[1, 3], mat[3, 3])
            yaw = atan2(mat[2, 1], mat[2, 2])
        else  # no singularity
            roll = atan2(-mat[2, 3], mat[3, 3])
            pitch = atan2(mat[1, 3], cos(roll) * mat[3, 3] - sin(roll) * mat[2, 3])
            yaw = atan2(-mat[1, 2], mat[1, 1])
        end
    else  # axis_order == :zyx
        singularity_present = all([abs(elem) < eps() for elem in [mat[1, 1] mat[2, 1]]])
        if singularity_present
            roll = 0
            pitch = atan2(-mat[3, 1], mat[1, 1])
            yaw = atan2(-mat[2, 3], mat[2, 2])
        else  # no singularity
            roll = atan2(mat[2, 1], mat[1, 1])
            sine_roll = sin(roll)
            cos_roll = cos(roll)
            pitch = atan2(-mat[3, 1], cos_roll * mat[1, 1] + sine_roll * mat[2, 1])
            yaw = atan2(sine_roll * mat[1, 3] - cos_roll * mat[2, 3],
                        cos_roll * mat[2, 2] - sine_roll * mat[1, 2])
        end
    end

    roll, pitch, yaw = convert_angle([roll, pitch, yaw], units)

    return [roll, pitch, yaw]
end

end # module transform
