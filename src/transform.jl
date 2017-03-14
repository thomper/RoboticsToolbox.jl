module transform

export rot2, rotx, roty, rotz
export trot2, trotx, troty, trotz
export se2, se3
export r2t, t2r
export rpy2r, rpy2t, rpy2jac, tr2rpy

typealias RealArray{T <: Real} Array{T,1}
typealias RealMatrix{T <: Real} Matrix{T}
typealias RealMatrix3D{T <: Real} Array{T,3}


function check_argument_units(units::Symbol)
    if !(units in (:rad, :deg))
        error("Expected :rad or :deg for units, got $units")
    end
end

function check_argument_axis_order(axis_order::Symbol)
    if !(axis_order in (:xyz, :zyx))
        error("Expected :xyz or :zyx for axis_order, got $axis_order")
    end
end

function convert_angle(θ::Union{Real, RealArray}, units::Symbol)
    check_argument_units(units)

    if units == :deg
            θ = deg2rad(θ)
    end

    return θ
end


# Rotation matrix generation
# -----------------------------------------------------------------------------

function rot2(θ::Real, units::Symbol=:rad)
    θ = convert_angle(θ, units)

    cos_θ = cos(θ)
    sin_θ = sin(θ)

    return [cos_θ -sin_θ;
            sin_θ cos_θ]
end

function rot_any(θ::Real, axis::Char, units::Symbol=:rad)
    θ = convert_angle(θ, units)

    cos_θ = cos(θ)
    sin_θ = sin(θ)

    if axis == 'x'
        return [1 0 0;
                0 cos_θ -sin_θ;
                0 sin_θ cos_θ]
    elseif axis == 'y'
        return [cos_θ 0 sin_θ;
                0 1 0;
                -sin_θ 0 cos_θ]
    elseif axis == 'z'
        return [cos_θ -sin_θ 0;
                sin_θ cos_θ 0;
                0 0 1]
    else
        error("Expected one of ('x', 'y', 'z') for argument axis, got $axis.")
    end
end

rotx(θ::Real, units::Symbol=:rad) = rot_any(θ, 'x', units)
roty(θ::Real, units::Symbol=:rad) = rot_any(θ, 'y', units)
rotz(θ::Real, units::Symbol=:rad) = rot_any(θ, 'z', units)


# Homogeneous transform generation
# -----------------------------------------------------------------------------

function r2t(rot_mat::Union{RealMatrix,RealMatrix3D})
    if ndims(rot_mat) == 3
        return cat(3, [r2t(rot_mat[:, :, i]) for i in 1:size(rot_mat, 3)]...)
    end

    if !(size(rot_mat) in ((2, 2), (3, 3)))
        println("ndims yo: $(ndims(rot_mat))")
        error("Expected array of size (2, 2) or (3, 3), instead had size $(size(rot_mat)).")
    end

    is_2d = size(rot_mat) == (2, 2)
    if is_2d
        return [rot_mat [0; 0]; [0 0 1]]
    else
        return [rot_mat [0; 0; 0]; [0 0 0 1]]
    end
end

function t2r(trans_mat::Union{RealMatrix, RealMatrix3D})
    if ndims(trans_mat) == 3
        return cat(3, [t2r(trans_mat[:, :, i]) for i in 1:size(trans_mat, 3)]...)
    end

    if !(size(trans_mat) in ((3, 3), (4, 4)))
        error("Expected array of size (3, 3) or (4, 4), instead had size $(size(rot_mat)).")
    end

    return trans_mat[1:end - 1, 1:end - 1]
end

trot2(θ::Real) = r2t(rot2(θ))

function trot_any(θ::Real, axis::Char, units::Symbol=:rad)
    rot_func = if axis == 'x'
                   rotx
               elseif axis == 'y'
                   roty
               elseif axis == 'z'
                   rotz
               else
                   error("Expected one of ('x', 'y', 'z') for axis, got $axis.")
               end
    return r2t(rot_func(θ, units))
end

trotx(θ::Real, units::Symbol=:rad) = trot_any(θ, 'x', units)
troty(θ::Real, units::Symbol=:rad) = trot_any(θ, 'y', units)
trotz(θ::Real, units::Symbol=:rad) = trot_any(θ, 'z', units)

function se2(x::Real, y::Real, θ::Real, units::Symbol=:rad)
    θ = convert_angle(θ, units)

    sin_θ = sin(θ)
    cos_θ = cos(θ)

    return [cos_θ -sin_θ x;
            sin_θ cos_θ y;
            0 0 1]
end

function se3(trans_mat::RealMatrix)
    if size(trans_mat) != (3, 3)
        error("Expected array of size (3, 3), instead had size $(size(trans_mat)).")
    end

    return [trans_mat[1:2, 1:2] [0; 0] trans_mat[1:2, 3];
            0 0 1 0;
            0 0 0 1]
end


# Conversion between roll/pitch/yaw and rotation matrices/homogeneous transforms
# ------------------------------------------------------------------------------

function rpy2r(roll::Real, pitch::Real, yaw::Real, units::Symbol=:rad, axis_order::Symbol=:xyz)
    check_argument_axis_order(axis_order)

    if axis_order == :xyz
        rot_func_a, rot_func_b, rot_func_c = (rotx, roty, rotz)
    elseif axis_order == :zyx
        rot_func_a, rot_func_b, rot_func_c = (rotz, roty, rotx)
    end

    return rot_func_a(roll, units) * rot_func_b(pitch, units) * rot_func_c(yaw, units)
end

rpy2t(roll::Real, pitch::Real, yaw::Real, units::Symbol=:rad, axis_order::Symbol=:xyz) = r2t(rpy2r(roll, pitch, yaw, units, axis_order))

function rpy2jac(roll::Real, pitch::Real, yaw::Real, units::Symbol=:rad)
    roll, pitch, yaw = convert_angle([roll, pitch, yaw], units)

    return [1 0 sin(pitch);
            0 cos(roll) (-cos(pitch) * sin(roll));
            0 sin(roll) (cos(pitch) * cos(roll))]
end

function tr2rpy(mat::Union{RealMatrix, RealMatrix3D}, units::Symbol=:rad, axis_order::Symbol=:xyz)
    if ndims(mat) == 3
        return cat(3, [tr2rpy(mat[:, :, i]) for i in 1:size(mat, 3)]...)
    end

    if !(size(mat) in ((3, 3), (4, 4)))
        error("Expected array of size (3, 3) or (4, 4), instead had size $(size(mat)).")
    end

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
