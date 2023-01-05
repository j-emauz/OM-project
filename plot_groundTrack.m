function [lon,lat]=plot_groundTrack(T_sat,n_orbits,theta_G0, y0,t0)
%PLOT_GROUNDTRACK Plots the ground track of a satellite given its initial
%conditions and the number of orbits it will perform
%
% The function uses the function ode_2bp to solve the 2-body problem
% using the ode113 algorithm. It then computes the longitude and
% latitude of the satellite using the position vector obtained from the
% integration. Finally, it plots the ground track of the satellite on
% top of a map of the Earth.
%
% INPUT
% - T_sat : period of the satellite [s]
% - n_orbits: number of orbits performed by the satellite
% - theta_G0: initial true anomaly of the satellite [rad]
% - y0 : initial conditions of the satellite (1x6 array)
% - t0 : initial time [s]
%
% OUTPUT
% - lon : longitude of the satellite at each time step [rad]
% - lat : latitude of the satellite at each time step [rad]
%
% k revolutions of the satellite
% m rotations of the planet
% Author:
% Name: Mariangela Testa, Oleksii Stepaniuk, João Emauz, Saverio Franzese
% Email: mariangela.testa@mail.polimi.it, oleksii.stepaniuk@mail.polimi.it,
% joao.emauz@mail.polimi.it, saverio.franzese@mail.polimi.it

mu_E = astroConstants(13);
omega_E=15.04 *pi/180/3600; %rad/s

tspan = linspace( 0, n_orbits*T_sat, 100000 );

options = odeset( 'RelTol', 2.3e-14, 'AbsTol', 1e-14 );
[ T, Y ] = ode113( @(t,y) ode_2bp(t,y,mu_E), tspan, y0, options );
    
delta=asin(Y(:,3)./norm(Y(:,1:3)));
alpha=atan2(Y(:,2), Y(:,1));
theta_G=theta_G0 + omega_E*(T-t0);

lon=wrapTo2Pi(alpha-theta_G);
lat=delta;

% graphics
figure;
img = imread('Earth_surface_plate.jpeg');
% image([0 0],[max(lon) max(lon)],img);
img= flip(img,1);  %# Flips the rows, making an upside-down image
hold on;
plot (lon.*(180/pi),lat.*(180/pi),'.',MarkerSize=0.7);
hold on
plot(lon(1)*180/pi,lat(1)*180/pi,'o',MarkerSize=9,LineWidth=2);
hold on
plot(lon(end)*180/pi,lat(end)*180/pi,'s',MarkerSize=9,LineWidth=2);
legend('Ground track','Start Position','End Position');
axis([-180 180 -90 90]);
xlabel('Longitude [deg]');
ylabel('Latitude [deg]');

end