function dy = perturbed_ode_2bp( ~, y, mu, J2, R)
%ode_2bp ODE system for the two-body problem (Keplerian motion)
%
% PROTOTYPE
% dy = perturbed_ode_2bp( t, y, mu, J2, R)
%
% INPUT:
% t[1] Time (can be omitted, as the system is autonomous) [T]
% y[6x1] State of the body ( rx, ry, rz, vx, vy, vz ) [ L, L/T ]
% mu[1] Gravitational parameter of the primary [L^3/T^2]
%
% OUTPUT:
% dy[6x1] Derivative of the state [ L/T^2, L/T^3 ]
%
% CONTRIBUTORS:
% Juan Luis Gonzalo Gomez
%
% VERSIONS
% 2018-09-26: First version
%
% -------------------------------------------------------------------------
% Position and velocity
r = y(1:3);
v = y(4:6);
% Distance from the primary
rnorm = norm(r);
% Set the derivatives of the state

a=3/2*J2*mu*R^2/rnorm^4*[r(1)/rnorm*(5*r(3)^2/rnorm^2-1); r(2)/rnorm*(5*r(3)^2/rnorm^2-1); r(3)/rnorm*(5*r(3)^2/rnorm^2-3)   ];
dy = [ v
(-mu/rnorm^3)*r + a];
end
