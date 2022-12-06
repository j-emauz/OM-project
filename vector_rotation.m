function v_rotated=vector_rotation(v,u,angle)
    v_rotated = v*cos(angle) + cross(u, v)*sin(angle) + u*(dot(u, v))*(1-cos(angle));
end 