function [dv_tot,VI,dt,TPAR] = dv_arc2(t1, t2, car_1, kep_2, mu)
% [dv_tot,t1,dt] = dv_calc(t1, t2, p1, p2, mu)
% This function calculates the total delta-v required for a spacecraft to go
% from one planet to another in a given time period.
% Inputs:
% t1: initial time (in days)
% t2: final time (in days)
% p1: initial planet
% p2: final planet
% mu: gravitational constant of central body
% Output:
% dv_tot: total delta-v required for the transfer

dt = t2 - t1;
dt = dt*24*3600;

[car_2, v_2] = par2car(kep_2(1), kep_2(2), kep_2(3), kep_2(4), kep_2(5), kep_2(6), mu);

[A,P,E,ERROR,VI,VF,TPAR,THETA] = lambertMR(car_1, car_2, dt, mu, 0, 0, 0 );

dv_tot = norm(v_2 - VF.');

end