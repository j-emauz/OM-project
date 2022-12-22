function [lon,lat]=groundTrack(T_sat,n_orbits,theta_G0, y0,t0)
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
mu_E = astroConstants(13);

tspan=linspace(0, n_orbits*T_sat, 100000);

options = odeset( 'RelTol', 2.3e-14, 'AbsTol', 1e-14 );
[ T, Y ] = ode113(@(t,y)ode_2bp(t,y,mu_E), tspan, y0, options );

delta=asin(Y(:,3)./vecnorm(Y(:,1:3),2,2));
alpha=atan2(Y(:,2), Y(:,1));
theta_G=theta_G0 + omega_E*(T-t0);


% for j=1:length(alpha)
%       lon(j)=alpha(j)-theta_G(j);
%     if alpha(j)-theta_G(j)>pi
%         lon(j)= alpha(j)-theta_G(j)-2*pi;
%     end 
%     if alpha(j)-theta_G(j)<-pi
%           lon(j)= alpha(j)-theta_G(j)+2*pi;  
%     end 
% end 

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
plot(lon_deg,lat_deg,'green', LineStyle='none',Marker='.', MarkerSize=7);
hold on;
plot(lon(1)*(180/pi),lat(1)*(180/pi),'o',MarkerSize=10,LineWidth=3, Color='green');
hold on;
plot(lon(end)*(180/pi),lat(end)*(180/pi),'s',MarkerSize=10,LineWidth=3, Color='red');
axis([-180 180 -90 90]);
xlabel('Longitude [deg]');
ylabel('Latitude [deg]');
legend('Ground track','Start','End');
end