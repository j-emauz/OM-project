% Set your initial time and state (Cartesian or Keplerian): t0, s0
% Set your integration time span: tspan
% Set physical parameters (e.g., gravitational parameter of the primary, J2, etc.)
% Set ODE solver options (e.g., RelTol, AbsTol): options
% Numerical integration of the equations of motion
clc
close all
clear

mu_E=astroConstants(13);
R_E=astroConstants(23);
J2=astroConstants(9);
kep0 = [7571, 0.01, 87.9 *pi/180, 180 *pi/180, 180 *pi/180, 0 *pi/180];
[car0, v] = par2car(kep0(1), kep0(2), kep0(3), kep0(4), kep0(5), kep0(6), mu_E);
T = 2*pi*sqrt( kep0(1)^3/mu_E ); % Orbital period [1/s]
tspan = linspace( 0, 100*T, 1000 );
options = odeset( 'RelTol', 1e-13, 'AbsTol', 1e-14 );
[ T, S ] = ode113( @(t,s) eq_motion( t, s, @(t,s) acc_pert_fun(t,s,mu_E,J2,R_E), mu_E,J2,R_E ), tspan, kep0, options );
% Analyse and plot the results
