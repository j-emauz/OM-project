function [dv_tot,VI,dt,TPAR] = dv_arc2(t1, t2, car_1, kep_2, mu)
% dv_arc2: Calculates the total delta-v required to transfer between two objects.
%
% INPUTS:
% t1 = Departure date [days]
% t2 = Arrival date [days]
% car_1 = Cartesian coordinates of departure position [km]
% kep_2 = Keplerian elements of arrival position
% mu = Gravitational parameter of the central body
%
% OUTPUTS:
% dv_tot = Total delta-v required for the transfer [km/s]
% VI = Initial velocity [km/s]
% dt = Time of flight [s]
% TPAR = Time of flight for minimum energy transfer [s]
%
% USAGE:
% [dv_tot,VI,dt,TPAR] = dv_arc2(t1, t2, car_1, kep_2, mu)
%
% Authors
% Name: Mariangela Testa, Oleksii Stepaniuk, João Emauz, Saverio Franzese
% Email: mariangela.testa@mail.polimi.it, oleksii.stepaniuk@mail.polimi.it,
% joao.emauz@mail.polimi.it, saverio.franzese@mail.polimi.it

% Compute time of flight
dt = t2 - t1;
dt = dt*24*3600;
% Convert Keplerian elements to Cartesian coordinates
[car_2, v_2] = par2car(kep_2(1), kep_2(2), kep_2(3), kep_2(4), kep_2(5), kep_2(6), mu);
% Compute minimum energy transfer
[A,P,E,ERROR,VI,VF,TPAR,THETA] = lambertMR(car_1, car_2, dt, mu, 0, 0, 0 );

% Calculate total delta-v
dv_tot = norm(v_2 - VF.');

end