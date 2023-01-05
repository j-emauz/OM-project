function v_rotated=vector_rotation(v,u,angle)
% vector_rotation: Rotate a Vector About an Axis
%
% This function rotates a vector "v" about an axis defined by a unit vector "u"
% by an angle "angle" (in radians). The output is the rotated vector "v_rotated".
%
% Usage: v_rotated = vector_rotation(v,u,angle)
%
% Inputs:
% v: vector to be rotated [km/s]
% u: unit vector defining rotation axis
% angle: rotation angle [rad]
%
% Outputs:
% v_rotated: rotated vector [km/s]
%
% Author:
% Name: Mariangela Testa, Oleksii Stepaniuk, João Emauz, Saverio Franzese
% Email: mariangela.testa@mail.polimi.it, oleksii.stepaniuk@mail.polimi.it,
% joao.emauz@mail.polimi.it, saverio.franzese@mail.polimi.it
    v_rotated = v*cos(angle) + cross(u, v)*sin(angle) + u*(dot(u, v))*(1-cos(angle));
end 