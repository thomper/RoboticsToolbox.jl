module transform

export rot2, rotx, roty, rotz
export trot2, trotx, troty, trotz
export se2, se3
export r2t, t2r
export rpy2r, rpy2t, rpy2jac, tr2rpy


# Rotation matrix generation
# -----------------------------------------------------------------------------

function rot2(θ::Number)
    cos_θ = cos(θ)
    sin_θ = sin(θ)

    return [cos_θ -sin_θ;
            sin_θ cos_θ]
end

function rot_any(θ::Number, axis::Char)
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

rotx = (θ::Number) -> rot_any(θ, 'x')
roty = (θ::Number) -> rot_any(θ, 'y')
rotz = (θ::Number) -> rot_any(θ, 'z')


# Homogeneous transform generation
# -----------------------------------------------------------------------------

# For some reason, using the Number parent class here gives an error for Array{Int64, 2}, e.g.:
# function r2t(rot_mat::Array{Number,2})
# Need to revisit this and add type restriction.  Likewise for all other functions in this
# module that accept arrays.
function r2t(rot_mat)
    if ~(size(rot_mat) in ((2, 2), (3, 3)))
        error("Expected array of size (2, 2) or (3, 3), instead had size $(size(rot_mat)).")
    end

    is_2d = size(rot_mat) == (2, 2)
    if is_2d
        return [rot_mat [0; 0]; [0 0 1]]
    else
        return [rot_mat [0; 0; 0]; [0 0 0 1]]
    end
end

# TODO Should be:
# function t2r(trans_mat::Array{Number,2})
function t2r(trans_mat)
    if ~(size(trans_mat) in ((3, 3), (4, 4)))
        error("Expected array of size (3, 3) or (4, 4), instead had size $(size(rot_mat)).")
    end

    return trans_mat[1:end - 1, 1:end - 1]
end

trot2 = (θ::Number) -> r2t(rot2(θ))

function trot_any(θ::Number, axis::Char)
    rot_func = if axis == 'x'
                   rotx
               elseif axis == 'y'
                   roty
               elseif axis == 'z'
                   rotz
               else
                   error("Expected one of ('x', 'y', 'z') for axis, got $axis.")
               end
    return r2t(rot_func(θ))
end

trotx = (θ::Number) -> trot_any(θ, 'x')
troty = (θ::Number) -> trot_any(θ, 'y')
trotz = (θ::Number) -> trot_any(θ, 'z')

function se2(x::Number, y::Number, θ::Number)
    sin_θ = sin(θ)
    cos_θ = cos(θ)

    return [cos_θ -sin_θ x;
            sin_θ cos_θ y;
            0 0 1]
end

# TODO Should be:
# function se3(trans_mat::Array{Number,2})
function se3(trans_mat)
    if size(trans_mat) != (3, 3)
        error("Expected array of size (3, 3), instead had size $(size(trans_mat)).")
    end

    return [trans_mat[1:2, 1:2] [0; 0] trans_mat[1:2, 3];
            0 0 1 0;
            0 0 0 1]
end


# Conversion between roll/pitch/yaw and rotation matrices/homogeneous transforms
# ------------------------------------------------------------------------------

rpy2r = (roll::Number, pitch::Number, yaw::Number) -> rotx(roll) * roty(pitch) * rotz(yaw)
rpy2t = (roll::Number, pitch::Number, yaw::Number) -> r2t(rpy2r(roll, pitch, yaw))

function rpy2jac(roll::Number, pitch::Number, yaw::Number)
    return [1 0 sin(pitch);
            0 cos(roll) (-cos(pitch) * sin(roll));
            0 sin(roll) (cos(pitch) * cos(roll))]
end

# TODO Should be:
# function tr2rpy(mat::Array{Number,2})
function tr2rpy(mat)
    if ~(size(mat) in ((3, 3), (4, 4)))
        error("Expected array of size (3, 3) or (4, 4), instead had size $(size(mat)).")
    end

    singularity_present = all([abs(elem) < eps() for elem in [mat[2, 3] mat[3, 3]]])
    if singularity_present
        roll = 0
        pitch = atan2(mat[1, 3], mat[3, 3])
        yaw = atan2(mat[2, 1], mat[2, 2])
    else
        roll = atan2(-mat[2, 3], mat[3, 3])
        pitch = atan2(mat[1, 3], cos(roll) * mat[3, 3] - sin(roll) * mat[2, 3])
        yaw = atan2(-mat[1, 2], mat[1, 1])
    end

    return roll, pitch, yaw
end

end # module transform
