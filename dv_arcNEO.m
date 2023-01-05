function [dv_tot,VF,v_2, car_2,t1,dt] = dv_arcNEO(t1, t2, p1, p2, mu)
% dv_arcNEO: Calculates the total delta-v required to transfer between planet and NEO.
%
% INPUTS:
% t1 = Departure date [days]
% t2 = Arrival date [days]
% p1, p2 = Departure and arrival objects
% mu = Gravitational parameter of the central body
%
% OUTPUTS:
% dv_tot = Total delta-v required for the transfer [km/s]
% VF = Final velocity [km/s]
% v_2 = Velocity of arrival planet at time t2 [km/s]
% car_2 = Cartesian coordinates of arrival object at time t2 [km]
% t1, dt = Modified departure date and time of flight [days]
%
% USAGE:
% [dv_tot,VF,v_2, car_2,t1,dt] = dv_arcNEO(t1, t2, p1, p2, mu)
%
% Authors
% Name: Mariangela Testa, Oleksii Stepaniuk, João Emauz, Saverio Franzese
% Email: mariangela.testa@mail.polimi.it, oleksii.stepaniuk@mail.polimi.it,
% joao.emauz@mail.polimi.it, saverio.franzese@mail.polimi.it

% Compute time of flight
dt = t2 - t1;
dt = dt*24*3600;
% Compute Keplerian elements at departure and arrival   
[kep_1,~] = uplanet(t1, p1);
[kep_2,~] = ephNEO(t2, p2);
% Convert Keplerian elements to Cartesian coordinates
[car_1, v_1] = par2car(kep_1(1), kep_1(2), kep_1(3), kep_1(4), kep_1(5), kep_1(6), mu);
[car_2, v_2] = par2car(kep_2(1), kep_2(2), kep_2(3), kep_2(4), kep_2(5), kep_2(6), mu);
% Compute delta-v required for the transfer using Lambert's problem
[A,P,E,ERROR,VI,VF,TPAR,THETA] = lambertMR(car_1, car_2, dt, mu, 0, 0, 0 );
% Calculate total delta-v
dv_tot =norm(VI.' - v_1);

end