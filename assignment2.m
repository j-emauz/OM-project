%% Assignment 2 - Planetary Explorer Mission
clc
clear 
close all

% Constants: 
AU=astroConstants(2);
mu_S=astroConstants(4);
mu_E=astroConstants(13);
R_E=astroConstants(23);
J2=astroConstants(9);
%% Orbit data and plot:

a=29032; %[km] 
e=0.4424; 
i=deg2rad(6.9341); %[rad]

%Assumption:
Om=90; %[deg]
om=70; %[deg]
theta0=0; %[deg], Starting at pericentre
kep0 = [a, e, i, Om *pi/180, om *pi/180, theta0 *pi/180];

%Conversion in Cartesian elements:
[r0, v0] = par2car(kep0(1), kep0(2), kep0(3), kep0(4), kep0(5), kep0(6), mu_E);

% Orbital period [sec]:
T = 2*pi*sqrt( kep0(1)^3/mu_E ); 

%radius of pericenter and apocenter:
rp=a*(1-e); %km 
ra=a*(1+e); %km

%Semi-Latus rectum:
p=a*(1-e^2); %km

%Specific angular momentum:
h=sqrt(p * mu_E); %kg*m^2/s


%Orbit integration:
y0 = [r0;v0];
tspan = linspace( 0, T,1000);
     % Set options for the ODE solver
options = odeset( 'RelTol', 1e-14, 'AbsTol', 1e-14 );
[ ~, Y ] = ode113( @(t,y) ode_2bp(t,y,mu_E), tspan, y0, options );
earth_sphere;
hold on;
grid on;
plot3( Y(:,1), Y(:,2), Y(:,3), 'r-','LineWidth',2);
title('Starting orbit');
legend('Earth Sphere', 'Starting orbit');

%% Ground Track of unperturbed 2BP:
omega_E=15.04 *pi/180/3600; %[rad/s]
theta_G0=0;
t0=0;

% 1 Orbit:
n_orbits=1;
[lon1,lat1]=groundTrack(T,n_orbits,theta_G0, y0,t0);
hold on; 
title('Ground Track unperturbed 2BP - 1 orbit');
% 1 day: 
n_orbits=24*3600/T;
[lon2,lat2]=groundTrack(T,n_orbits,theta_G0, y0,t0);
hold on; 
title('Ground Track unperturbed 2BP - 1 day');

% 10 days: 
n_days=10;
n_orbits=n_days*24*3600/T;
[lon3,lat3]=groundTrack(T,n_orbits,theta_G0, y0,t0);
hold on; 
title('Ground Track unperturbed 2BP - 10 days');
hold off;
%% Repeated Ground Track:
% k = revolutions of the satellite;
% m = rotations of the planet;

%Using the ratio k/m assigned: 
k=7;
m=4;
ratio=k/m;
n_orbits=k;
n_SC=omega_E*ratio; %[rad/s]
a_repeating=(mu_E/n_SC^2)^(1/3);
T_repeating=2*pi/omega_E*1/ratio;

[lon4,lat4]=groundTrack(T_repeating,n_orbits,theta_G0, y0,t0);
hold on; 
title('Repeated Ground Track unperturbed 2BP - Given ratio k/m=7/4');
hold off;

% % 1 orbit: 
% k=1; 
% [a_repeating_1, k_1, m_1]=RepeatingGroundTracks(k,omega_E,T,mu_E,1);
% T_repeating_1=2*pi*sqrt(a_repeating_1^3/mu_E);
% [lon1,lat1]=groundTrack(T_repeating_1,k_1,theta_G0, y0,t0);
% hold on; 
% title('Repeated Ground Track unperturbed 2BP - 1 orbit');
% hold off;
% % 1 day: 
% m=1;
% [a_repeating_2, k_2,m_2]=RepeatingGroundTracks(m,omega_E,T,mu_E,0);
% T_repeating_2=2*pi*sqrt(a_repeating_2^3/mu_E);
% [lon2,lat2]=groundTrack(T_repeating_2,k_2,theta_G0, y0,t0);
% hold on; 
% title('Repeated Ground Track unperturbed 2BP - 1 day');
% hold off;
% % 10 days:
% m=10;
% [a_repeating_3, k_3, m_3]=RepeatingGroundTracks(m,omega_E,T,mu_E,0);
% T_repeating_3=2*pi*sqrt(a_repeating_3^3/mu_E);
% [lon3,lat3]=groundTrack(T_repeating_3,k_3,theta_G0, y0,t0);
% hold on; 
% title('Repeated Ground Track unperturbed 2BP - 10 days');
% hold off;


%% Perturbed 2BP:

% Data:
AMR= 10.000; %[m^2/kg]
P_sun=4.5*10^-6; %[N/m^2]
Cr=1; 

% Earth orbit around the Sun:
initial_date=[2022,03,21,12,0,0]; %Spring Equinox
initial_time=date2mjd2000(initial_date);
[kep_E,~]=uplanet(initial_time,3);
[r0,v0] = par2car(kep_E(1),kep_E(2),kep_E(3),kep_E(4),kep_E(5),kep_E(6),mu_S);
y0 = [r0;v0];
% T=2*pi*sqrt(kep_E(1)^3/mu_S);
% tspan = linspace( 0, T,1000);

% Set options for the ODE solver
options = odeset( 'RelTol', 1e-14, 'AbsTol', 1e-14 );
[ t, Y_earth ] = ode113( @(t,y) ode_2bp(t,y,mu_S), tspan, y0, options );

r_E_S=Y_earth(:,1:3);
%r_sc_E=Y(:,1:3);
%r_sc_Sun_ecliptic=r_E_S;%r_E_S + r_sc_E
% r_sc_Sun=zeros(length(tspan),3);
    
%tilt = deg2rad(23.45);
%R1_eps = [1   0   0;                %ECI -> Sun-Centered-Ecliptic
         % 0  cos(tilt)   sin(tilt);
          %0  -sin(tilt)  cos(tilt)];
%r_sc_Sun_ECI= R1_eps'*(r_sc_Sun_ecliptic); 

%for j=1:size(Y(:,1))
%    [~,~,~,~,~,theta_sc,~]=car2par(Y(j,1:3)',Y(j,4:6)',mu_E);
%     R1_i=[1 0 0; 0 cos(i) sin(i); 0 -sin(i) cos(i)];
%     R3_Om=[cos(Om) sin(Om) 0; -sin(Om) cos(Om) 0; 0 0 1];
%     R3_om_th=[cos(om+theta_sc) sin(om+theta_sc) 0; -sin(om+theta_sc) cos(om+theta_sc) 0; 0 0 1];
    
    %r_sc_Sun_ECI(j,:)= R1_eps*(r_sc_Sun_ecliptic(j,:)');  
%end 


%Orbit propagation with Cartesian Coordinates:
tspan_pert = linspace( 0, 100*T, 1000 );
options = odeset( 'RelTol', 1e-13, 'AbsTol', 1e-14 );
[ T_J2, S_J2 ] = ode113(@(t,s) perturbed_ode_2bp_SRP(t,s, mu_E, J2, R_E), tspan_pert, y0, options );

for k=1:size(S_J2,1)
   [a_p(k),e_p(k),i_p(k),OM_p(k),om_p(k),th_p(k),~]=car2par(S_J2(k,1:3)', S_J2(k,4:6)',mu_E);
%[a,e,i,OM,om,th, ee]=car2par(rr,vv,mu)
%@(t_tt, y_tt) ode_2bp(t_tt, y_tt, mu), tspan, y0_t, options);
end

%Orbit propagation in Keplerian Elements using Gauss' planetary equations:
tspan = linspace( 0, 100*T, 1000 );
options = odeset( 'RelTol', 1e-13, 'AbsTol', 1e-14);
[ T_Gauss, S_Gauss ] = ode113(@(t,s) eq_motion(t,s, @(t,s) acc_pert_fun_J2_SRP(t,s,mu_E,J2, R_E), mu_E), tspan_pert, kep0, options );

% Semi-Major Axis, a
N = 100;
amean = movmean(S_Gauss(:, 1), N);
hold on;
plot(T_Gauss/T, amean);
plot(T_J2/T,a_p);
plot(T_Gauss/T,S_Gauss(:,1));
legend('average','Cartesian','Gauss')

a_error = (abs(a_p-S_Gauss(:,1)'))/kep0(1);

figure
grid on
semilogy(T_Gauss/T,a_error);

% Eccentricity, e
figure;
N = 500;
emean = movmean(S_Gauss(:, 2), N);
plot(T_Gauss/T,S_Gauss(:,2));
hold on
plot(T_J2/T,e_p);
plot(T_Gauss/T, emean);
legend('Gauss','Cartesian','average')

e_error = (abs(S_Gauss(:,2)'-e_p));
figure
grid on
semilogy(T_Gauss/T,e_error);

% Inclination, i
figure;
N = 500;
imean = movmean(S_Gauss(:, 3), N);
plot(T_Gauss/T,S_Gauss(:,3));
hold on
plot(T_J2/T,i_p);
plot(T_Gauss/T, imean);
legend('Gauss','Cartesian','average')

i_error = (abs(S_Gauss(:,3)'-i_p))/(2*pi);
figure
grid on
semilogy(T_Gauss/T,i_error);

% Right Ascension of the Ascending Node, OM
figure;
N = 500;
OMmean = movmean(S_Gauss(:, 4), N);
plot(T_Gauss/T,S_Gauss(:,4));
hold on
plot(T_J2/T,OM_p);
plot(T_Gauss/T, OMmean);
legend('Gauss','Cartesian','average')

OM_error = (abs(S_Gauss(:,4)'-OM_p))/(2*pi);
figure
grid on
semilogy(T_Gauss/T,OM_error);

% Argument of pericenter, om
figure;
N = 500;
ommean = movmean(S_Gauss(:, 5), N);
plot(T_Gauss/T,S_Gauss(:,5));
hold on
plot(T_J2/T,om_p);
plot(T_Gauss/T, ommean);
legend('Gauss','Cartesian','average')

om_error = (abs(S_Gauss(:,5)'-om_p))/(2*pi);
figure
grid on
semilogy(T_Gauss/T,om_error);

% True Anomaly, theta
figure;
S_Gauss(:, 6) = wrapTo2Pi(S_Gauss(:,6));
N = 500;
thmean = movmean(S_Gauss(:, 6), N);
plot(T_Gauss/T,S_Gauss(:,6));
hold on
plot(T_J2/T,th_p);
plot(T_Gauss/T, thmean);
legend('Gauss','Cartesian','average')

th_error = (abs(S_Gauss(:,6)'-th_p))/(2*pi);
figure
grid on
semilogy(T_Gauss/T,th_error);
