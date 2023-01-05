function acc_pert_vec = acc_pert_fun( t, s, mu,J2,R)
% acc_pert_fun: Perturbing Acceleration in a Given RSW Frame
%
% This function evaluates the perturbing acceleration of an object in a given
% rotating, synodic, and wavering (RSW) frame. The input parameters are the
% time "t", state vector "s", gravitational constant "mu", J2 coefficient, and
% radius "R" of the attractor body. The output parameter is the perturbing
% acceleration vector "acc_pert_vec".
%
% The state vector "s" consists of the following elements:
% s(1) = semi-major axis [km]
% s(2) = eccentricity
% s(3) = inclination [rad]
% s(4) = right ascension of the ascending node [rad]
% s(5) = argument of perigee [rad]
% s(6) = true anomaly [rad]
%
% Usage: acc_pert_vec = acc_pert_fun( t, s, mu,J2,R)
%
% INPUT:
% t: time [s]
% s: state vector [km, -, rad, rad, rad, rad]
% mu: gravitational constant [km^3/s^2]
% J2: coefficient
% R: radius of attractor body [km]
%
% OUTPUT:
% acc_pert_vec: perturbing acceleration vector [km/s^2]
%
    p = s(1)*(1-s(2)^2);
    r = p/(1 + s(2)*cos(s(6)));
    v = sqrt(2*mu/r - mu/s(1));
    h = sqrt(p*mu);
    
    acc_pert_vec = -3/2*J2*mu*R^2/r^4.*[1-3*(sin(s(3)))^2*sin((s(6)+s(5)))^2; (sin(s(3)))^2*sin(2*(s(6)+s(5))); sin(2*s(3))*sin((s(6)+s(5)))];
    
end