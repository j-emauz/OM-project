function [lon,lat]=groundTrack_perturbed(T_sat,n_orbits,theta_G0,initial_date, AMR, Cr,kep0,t0)
% INPUTS
% r, at the initial time (either in Cartesian or Keplerian elements)
% theta_G0, Longitude of Greenwich meridian at initial time
% t, Vector of times at which the ground track will be computed


% OUTPUTS
% alpha: right ascension in Earth Centred Equatorial Inertial frame
% delta: declination in Earth Centred Equatorial Inertial frame
% lon: longitude with respect to rotating Earth (0 deg at Greenwich meridian)
% lat: latitude with respect to rotating Earth

omega_E=15.04 *pi/180/3600; %rad/s
R_E=astroConstants(23);
mu_E = astroConstants(13);
J2=astroConstants(9);

tspan_pert=linspace(0, n_orbits*T_sat, 100000);

options = odeset( 'RelTol', 2.3e-14, 'AbsTol', 1e-14 );
% [ T, Y ] = ode113(@(t,y)ode_2bp(t,y,mu_E), tspan, y0, options );
[ T_Gauss, S_Gauss ] = ode113(@(t,s) eq_motion(t,s, @(t,s) acc_pert_fun_J2_SRP(t,s,mu_E,J2, R_E, initial_date, AMR, Cr, 2), mu_E), tspan_pert, kep0, options );

for j=1:size(S_Gauss,1)
    [r_Gauss(j,:), ~] = par2car(S_Gauss(j,1), S_Gauss(j,2), S_Gauss(j,3), S_Gauss(j,4), S_Gauss(j,5), S_Gauss(j,6), mu_E);
end

delta=asin(r_Gauss(:,3)./vecnorm(r_Gauss(:,1:3),2,2));
alpha=atan2(r_Gauss(:,2), r_Gauss(:,1));
theta_G=theta_G0 + omega_E*(T_Gauss-t0);


lon=wrapToPi(alpha-theta_G);
lat=delta;
lon_deg=lon.*(180/pi);
lat_deg=lat.*(180/pi);

%Graphics:
figure;
hold on;
img = imread('land_shallow_topo_2048.tif');
% img= flip(img,1);  % Flips the rows, making an upside-down image;
image([-180 180],[-90 90],flip(img));
hold on;
plot(lon_deg,lat_deg,'green', LineStyle='none',Marker='.', MarkerSize=4);
hold on;
plot(lon(1)*(180/pi),lat(1)*(180/pi),'o',MarkerSize=10,LineWidth=3, Color='green');
hold on;
plot(lon(end)*(180/pi),lat(end)*(180/pi),'s',MarkerSize=10,LineWidth=3, Color='red');
axis([-180 180 -90 90]);
xlabel('Longitude [deg]');
ylabel('Latitude [deg]');
legend('Perturbed ground track','Start','End');
end

