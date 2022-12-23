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
Om=35; %[deg]
om=25; %[deg]
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
legend('', 'Starting orbit');




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


%% Perturbed 2BP:

% Data:
AMR= 10.000; %[m^2/kg]
P_sun=4.5*10^-6; %[N/m^2]
Cr=1; 
initial_date=[2022,03,21,12,0,0];


% Orbit propagation with Cartesian Coordinates:

tspan_pert = linspace( 0, 1000*T, 100000 );
options = odeset( 'RelTol', 1e-13, 'AbsTol', 1e-14 );

% computational time
t = tic;
tCPU_Cartesian_start = cputime;
[ T_J2, S_perturbed ] = ode113(@(t,s) perturbed_ode_2bp_SRP(t,s, mu_E, J2, R_E, initial_date, AMR, Cr, 2), tspan_pert, y0, options );
tCPU_Cartesian = cputime - tCPU_Cartesian_start;
time_computation = toc (t);

% plotting perturbed orbit (1)
figure
earth_sphere
hold on
plot3( S_perturbed(:,1), S_perturbed(:,2), S_perturbed(:,3), 'r-','LineWidth',1);
hold off

% plotting perturbed orbit (2)
% view([270 45])
figure
comet3(S_perturbed(:,1), S_perturbed(:,2), S_perturbed(:,3))
title('Perturbed orbit using Cartesian coordinates');
hold off

% converting Cartesian coordinates into keplerian parameters for comparison
% in the following plots
for k=1:size(S_perturbed,1)
   [a_p(k),e_p(k),i_p(k),OM_p(k),om_p(k),th_p(k),~]=car2par(S_perturbed(k,1:3)', S_perturbed(k,4:6)',mu_E);
%[a,e,i,OM,om,th, ee]=car2par(rr,vv,mu)
%@(t_tt, y_tt) ode_2bp(t_tt, y_tt, mu), tspan, y0_t, options);
end

%% Orbit propagation in Keplerian Elements using Gauss' planetary equations:
options = odeset( 'RelTol', 1e-13, 'AbsTol', 1e-14);

% computational time
tCPU_Gauss_start = cputime;
[ T_Gauss, S_Gauss ] = ode113(@(t,s) eq_motion(t,s, @(t,s) acc_pert_fun_J2_SRP(t,s,mu_E,J2, R_E, initial_date, AMR, Cr, 2), mu_E), tspan_pert, kep0, options );
tCPU_Gauss = cputime - tCPU_Gauss_start;


for j=1:size(S_Gauss,1)
    [r_Gauss(j,:), ~] = par2car(S_Gauss(j,1), S_Gauss(j,2), S_Gauss(j,3), S_Gauss(j,4), S_Gauss(j,5), S_Gauss(j,6), mu_E);
end

% hold on
% plot3( r_Gauss(:,1), r_Gauss(:,2), r_Gauss(:,3));
% hold off
figure
 earth_sphere
 x = r_Gauss(:,1);
 y = r_Gauss(:,2); 
 z = r_Gauss(:,3);
view(3);
hold on; grid on;
c = colorbar;
c.Label.String = 'Time [s]';
 color = jet(length(T_Gauss)-1);
 for k = 2:size(T_Gauss,1)
     plot3(x(k-1:k),y(k-1:k),z(k-1:k),'color',color(k-1,:))
%     axis([0 1800.7 -1.7 1.7 -1.3 1.3])
    
 end

% plotting perturbed orbit (2)
% view([270 45])
% comet3(r_Gauss(:,1), r_Gauss(:,2), r_Gauss(:,3))
% title('Perturbed orbit using Cartesian coordinates');
% hold off

% Semi-Major Axis, a
figure
N = 100;
amean = movmean(S_Gauss(:, 1), N);
hold on;
plot(T_J2/T,a_p);
plot(T_Gauss/T,S_Gauss(:,1));
plot(T_Gauss/T, amean);
legend('Cartesian','Gauss','average');
title('semi-major axis');
a_error = (abs(a_p-S_Gauss(:,1)'))/kep0(1);

figure
grid on
semilogy(T_Gauss/T,a_error);
title('Semi-major axis error');

% Eccentricity, e
figure;
N = 100;
emean = movmean(S_Gauss(:, 2), N);
emean2 = movmean(S_Gauss(:,2),1000000);
plot(T_Gauss/T,S_Gauss(:,2));
title('eccentricity');
hold on
plot(T_J2/T,e_p);
plot(T_Gauss/T, emean);
plot(T_Gauss/T, emean2)
legend('Gauss','Cartesian','short-term average', 'secular average')

e_error = (abs(S_Gauss(:,2)'-e_p));
figure
grid on
semilogy(T_Gauss/T,e_error);
title('eccentricity error');

% Inclination, i
figure;
N = 100;
imean = movmean(S_Gauss(:, 3), N);
plot(T_Gauss/T,S_Gauss(:,3));
hold on
plot(T_J2/T,i_p);
plot(T_Gauss/T, imean);
title('inclination');
legend('Gauss','Cartesian','average')

i_error = (abs(S_Gauss(:,3)'-i_p))/(2*pi);
figure
grid on
semilogy(T_Gauss/T,i_error);
title('inclination error');

% Right Ascension of the Ascending Node, OM
figure;
N = 100;
OMmean = movmean(S_Gauss(:, 4), N);
plot(T_Gauss/T,S_Gauss(:,4));
hold on
plot(T_J2/T,OM_p);
plot(T_Gauss/T, OMmean);
legend('Gauss','Cartesian','average')
title('Right Ascension of the Ascending Node');

OM_error = (abs(S_Gauss(:,4)'-OM_p))/(2*pi);
figure
grid on
semilogy(T_Gauss/T,OM_error);
title('Right Ascension of the Ascending Node error');

% Argument of pericenter, om
figure;
N = 100;
ommean = movmean(S_Gauss(:, 5), N);
plot(T_Gauss/T,S_Gauss(:,5));
hold on
plot(T_J2/T,om_p);
plot(T_Gauss/T, ommean);
legend('Gauss','Cartesian','average')
title('argument of pericenter');

om_error = (abs(S_Gauss(:,5)'-om_p))/(2*pi);
figure
grid on
semilogy(T_Gauss/T,om_error);
title('argument of pericenter error');

% True Anomaly, theta
figure;
S_Gauss(:, 6) = wrapTo2Pi(S_Gauss(:,6));
N = 100;
thmean = movmean(S_Gauss(:, 6), N);
plot(T_Gauss/T,S_Gauss(:,6));
hold on
plot(T_J2/T,th_p);
plot(T_Gauss/T, thmean);
legend('Gauss','Cartesian','average')
title('true anomaly');

th_error = (abs(S_Gauss(:,6)'-th_p))/(2*pi);
figure
grid on
semilogy(T_Gauss/T,th_error);
title('true anomaly error');

%% point 7 - Comparison with real data

sc_ephemeris=load("056B_3y_matrix.mat");

%e=sc_ephemeris.B2hours.EC;
% a=sc_ephemeris.B2hours.A;
% i=sc_ephemeris.B2hours.IN.*pi/180;
% OM=wrapToPi(sc_ephemeris.B2hours.OM.*pi/180);
% om=wrapToPi(sc_ephemeris.B2hours.W.*pi/180);
% theta=wrapToPi(sc_ephemeris.B2hours.TA.*pi/180);
% date=sc_ephemeris.B2hours.CalendarDateTDB; %date

e=sc_ephemeris.B3y.EC;
a=sc_ephemeris.B3y.A;
i=sc_ephemeris.B3y.IN.*pi/180;
OM=wrapToPi(sc_ephemeris.B3y.OM.*pi/180);
om=wrapToPi(sc_ephemeris.B3y.W.*pi/180);
theta=wrapToPi(sc_ephemeris.B3y.TA.*pi/180);
date=sc_ephemeris.B3y.CalendarDateTDB; %date
date0=date(1);

kep0_sc=[a(1),e(1),i(1),OM(1),om(1),theta(1)];

date_0=[2020,12,22,0,0,0];
t0=date2mjd2000(date_0);
date_final=[2023,1,19,0,0,0];
t_final=date2mjd2000(date_final);

%t_span_sc=t0*24*3600:2*3600:t_final*24*3600;
t_span_sc=t0*24*3600:24*3600:t_final*24*3600;

for j=1:length(e)
    [rr_sc_ephemeris(:,j), vv_sc_ephemeris(:,j)] = par2car(a(j), e(j), i(j), OM(j), om(j), theta(j), mu_E);
end

earth_sphere
hold on
plot3( rr_sc_ephemeris(1,:), rr_sc_ephemeris(2,:), rr_sc_ephemeris(3,:),'r-','LineWidth',1);
hold off

%% 7.c
close all
AMR = 0.1;
Cr = 1;
T_sc = 2*pi*sqrt( kep0_sc(1)^3/mu_E )/3600; 
% tspan_pert = linspace( 0, 100*T_sc, 10000 );
[ T_sc_Gauss, S_sc_Gauss ] = ode113(@(t,s) eq_motion(t,s, @(t,s) acc_pert_fun_J2_SRP(t,s,mu_E,J2, R_E, date_0, AMR, Cr, 2), mu_E), t_span_sc, kep0_sc, options );
%[rr, vv] = par2car(kep0_sc(1), kep0_sc(2), kep0_sc(3), kep0_sc(4), kep0_sc(5), kep0_sc(6), mu_E);
%kep0_sc = [rr; vv];
%[ T_sc_Gauss, S_sc_Gauss ] = ode113(@(t,s) perturbed_ode_2bp_SRP(t,s, mu_E,J2, R_E, date_0, AMR, Cr, 2), t_span_sc, kep0_sc, options );
% Semi-Major Axis, a
% N = 100;
% amean = movmean(S_Gauss(:, 1), N);
% plot(T_sc_Gauss/T_sc, amean);
%% a
figure;
hold on;
plot(T_sc_Gauss/T_sc,S_sc_Gauss(:,1));
plot(T_sc_Gauss/T_sc,a);

legend('Gauss','Ephemerides');
title('semi-major axis');
hold off


%% Eccentricity, e
figure;
% N = 500;
% emean = movmean(S_Gauss(:, 2), N);
plot(T_sc_Gauss/T_sc,S_sc_Gauss(:,2));
title('eccentricity');
hold on
plot(T_sc_Gauss/T_sc,e);
% plot(T_Gauss/T, emean);
legend('Gauss','Ephemeris')

%% Inclination, i
figure;
% N = 500;
% imean = movmean(S_Gauss(:, 3), N);
plot(T_sc_Gauss/T_sc,S_sc_Gauss(:,3));
hold on
plot(T_sc_Gauss/T_sc,i);
% plot(T_Gauss/T, imean);
title('inclination');
legend('Gauss','Ephemeris')

%% Right Ascension of the Ascending Node, OM
figure;
% N = 500;
% OMmean = movmean(S_Gauss(:, 4), N);
plot(T_sc_Gauss/T_sc,S_sc_Gauss(:,4));
hold on
plot(T_sc_Gauss/T_sc,OM);
% plot(T_Gauss/T, OMmean);
legend('Gauss','Ephemeris')
title('Right Ascension of the Ascending Node');

%% Argument of pericenter, om
figure;
% N = 500;
% ommean = movmean(S_Gauss(:, 5), N);
plot(T_sc_Gauss/T_sc,S_sc_Gauss(:,5));
hold on
plot(T_sc_Gauss/T_sc,om);
% plot(T_sc_Gauss/T_sc, ommean);
legend('Gauss','Ephemeris')
title('argument of pericenter');



%% True Anomaly, theta
figure;
S_sc_Gauss(:, 6) = wrapTo2Pi(S_sc_Gauss(:,6));
% N = 500;
% thmean = movmean(S_Gauss(:, 6), N);
plot(T_sc_Gauss/T_sc,S_sc_Gauss(:,6));
hold on
plot(T_sc_Gauss/T_sc,theta);
% plot(T_sc_Gauss/T, thmean);
legend('Gauss','Ephemeris')
title('true anomaly');
