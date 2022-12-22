function dy = perturbed_ode_2bp( t, y, mu, J2, R)
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

Cr = 1;
AU = astroConstants(2);
AMR = 10;
mu_S = astroConstants(4);
P_sun = 4.5*10^-6;

 initial_date=[2022,03,21,12,0,0];
 initial_time=date2mjd2000(initial_date);
 t_days=t/3600/24/365.25; %days since initial time
 [kep_E,~]=uplanet(initial_time+t_days,3);
 [r_E_S,~] = par2car(kep_E(1),kep_E(2),kep_E(3),kep_E(4),kep_E(5),kep_E(6),mu_S);

    
  tilt = deg2rad(23.45);
  R_eci2sce = [1   0   0;                %ECI -> Sun-Centered-Ecliptic
        0  cos(tilt)   sin(tilt);
        0  -sin(tilt)  cos(tilt)];
  r_sc_Sun= R_eci2sce'*(r_E_S); 
  acc_pert_SRP= P_sun*(AU^2)/(norm(r_sc_Sun)^2)*Cr*AMR;
  acc_pert_vec_SRP=-acc_pert_SRP*(r_sc_Sun./norm(r_sc_Sun));

a=3/2*J2*mu*R^2/rnorm^4*[r(1)/rnorm*(5*r(3)^2/rnorm^2-1); r(2)/rnorm*(5*r(3)^2/rnorm^2-1); r(3)/rnorm*(5*r(3)^2/rnorm^2-3)   ];
dy = [ v
(-mu/rnorm^3)*r + a + acc_pert_vec_SRP];
end
