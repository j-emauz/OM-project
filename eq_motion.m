function ds = eq_motion( t, s, acc_pert_fun, parameters )
% Evaluate the perturbing accelerations
acc_pert_vec = acc_pert_fun( t, s );
% Evaluate the equations of motion (Cartesian or Keplerian),
% function of t, s, acc_pert_vec, and parameters
ds = ...;

r = y(1:3);
v = y(4:6);
% Distance from the primary
rnorm = norm(r);
% Set the derivatives of the state
dy = [ v
(-mu/rnorm^3)*r ];




end