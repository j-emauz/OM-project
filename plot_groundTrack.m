function [lon,lat]=plot_groundTrack(T_sat,n_orbits,theta_G0, y0,t0)

% k revolutions of the satellite
% m rotations of the planet

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