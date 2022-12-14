clear all
clc
close all

% constants
mu_S = astroConstants(4);
mu_E= astroConstants(13);
AU= astroConstants(2);

% dates
dept1=[2028,01,30,0,0,0]; %GIVEN
arrt2=[2062,07,28,0,0,0]; % GIVEN latest arrival on the asteroid

% synodic period
a_Saturn=9.5826*AU;
T_Saturn = 2*pi*sqrt( a_Saturn^3/mu_S ); % Orbital period [1/s]
a_Earth=149598023;
T_Earth=2*pi*sqrt( a_Earth^3/mu_S );
Tsyn_ES=T_Saturn*T_Earth/(abs(T_Saturn-T_Earth));

dept2_mjd2000 = date2mjd2000(dept1)+5*Tsyn_ES/3600/24;
dept2=mjd20002date(dept2_mjd2000);
% dept2=[2030,01,30,0,0,0];  % random

arrt1=[2043,07,28,0,0,0]; % random value


GA_time_min=[2032,01,30,0,0,0]; %assumption for fly-by window
GA_time_max=[2033,01,30,0,0,0]; %assumption for fly-by window

%date2mjd2000
%conversion to Modern Julian Date 2000
t11 = date2mjd2000(dept1);
t12 = date2mjd2000(dept2);
t21 = date2mjd2000(arrt1);
t22 = date2mjd2000(arrt2);
tGA1= date2mjd2000(GA_time_min);
tGA2= date2mjd2000(GA_time_max);

tspan1 = linspace(t11, t12, 100); %departure window
tspan2 = linspace(t21, t22, 100); %arrival window
% was 500
tspanGA=linspace(tGA1, tGA2, 100); % fly-by window

p1=3; %Earth
p2=6; %Saturn

for i=1:length(tspan1)
    for j=1:length(tspanGA)
        [dv_1(i,j),t1(i,j),dt] = dv_calc(tspan1(i), tspanGA(j), p1, p2, mu_S);
% velocity to escape from the parking orbit
        for k=1:length(tspan2)
            
        [dv_2(i,j,k),t2(i,j,k),dt] = dv_calc(tspanGA(i), tspan2(k), p1, p2, mu_S);
        dv_tot(i,j,k)=dv_2(i,j,k)+dv_1(i,j);
        end
    end
end

m1=min(dv_tot);
m2=min(min(dv_tot));
m3=min(min(min(dv_tot)));

% [x,y,z]=find(dv_tot==m3);
[x,y,z] = ind2sub(size(dv_tot),find(dv_tot==m3));
DepartureTime_min_dv=mjd20002date(tspan1(x));
FlyBy_Time_min_dv=mjd20002date(tspanGA(y));
ArrivalTime_min_dv=mjd20002date(tspan2(z));


%% OPTIMAL SOLUTION - 1st ARC
[kepEarth_dep,~] = uplanet(tspan1(x), p1); % position of Earth at departure time
[kepEarth_fb,~] = uplanet(tspanGA(y), p1); % position of Earth at fly by time
[kepSaturn_dep,~] = uplanet(tspan1(x), p2); % position of Saturn at departure time
[kepSaturn_fb,~]= uplanet(tspanGA(y), p2); % position of Saturn at fly by
[kepSaturn_arr,~]= uplanet(tspan2(z), p2); % position of Saturn at arrival

[kepNEO_fb,~,~] = ephNEO(tspanGA(y),86);
[kepNEO_arr,mass,M] = ephNEO(tspan2(z),86);

[car_Earth_dep, v_Earth_dep] = par2car(kepEarth_dep(1), kepEarth_dep(2), kepEarth_dep(3), kepEarth_dep(4), kepEarth_dep(5), kepEarth_dep(6), mu_S);
[car_Earth_fb, v_Earth_fb] = par2car(kepEarth_fb(1), kepEarth_fb(2),kepEarth_fb(3), kepEarth_fb(4), kepEarth_fb(5), kepEarth_fb(6), mu_S);
[car_Saturn_dep, v_Saturn_dep] = par2car(kepSaturn_dep(1),kepSaturn_dep(2),kepSaturn_dep(3),kepSaturn_dep(4), kepSaturn_dep(5), kepSaturn_dep(6), mu_S);
[car_Saturn_fb, v_Saturn_fb] = par2car(kepSaturn_fb(1),kepSaturn_fb(2),kepSaturn_fb(3),kepSaturn_fb(4), kepSaturn_fb(5), kepSaturn_fb(6), mu_S);

% car is used into the Lambert function and in y0

% position of planets at BEGINNING of departure and arrival windows

[kep_21,~] = uplanet(tGA1, p2);

[car_21, v_21] = par2car(kep_21(1), kep_21(2), kep_21(3), kep_21(4), kep_21(5), kep_21(6), mu_S);

% calculating transfer orbit BETWEEN THE OPTIMAL POINTS
[A,P,E,ERROR,VI,VF,TPAR,THETA] = lambertMR(car_Earth_dep, car_Saturn_fb, tspanGA(y)-tspan1(x), mu_S, 0, 2, 2 );
%[A,P,E,ERROR,VI,VF,TPAR,THETA] = lambertMR(RI,RF,TOF,MU,orbitType,Nrev,Ncase,optionsLMR)
%                 [kepNEO,mass,M,d] = ephNEO(time,86);
%                 [kepEarth,~] = uplanet(mjd2000, 3);
%                 [kepSaturn,~] = uplanet(mjd2000, 6);



%% PROPAGATION
options = odeset( 'RelTol', 1e-13, 'AbsTol', 1e-14 );

% 1 Tranfer orbit arc
tspan = linspace(0, tspanGA(y)-tspan1(x), 1000); %ToF
y0_t = [car_Earth_dep; VI.'];
% velocity of the transfer orbit in the first optimal point
[ t_t1, y_t1 ] = ode45(@(t_t, y_t) ode_2bp(t_t, y_t, mu_S), tspan, y0_t, options);

% 2 Earth orbit during transfer
tspan = linspace(0, tspanGA(y)-tspan1(x), 1000);
y0_1 = [car_Earth_dep; v_Earth_dep]; % optimal point on Earth orbit
[ t_1, y_1 ] = ode45(@(t_1, y_1) ode_2bp(t_1, y_1, mu_S), tspan, y0_1, options);

% 3 Earth orbit
a = kepEarth_dep(1); % Semi-major axis [km]
T = 2*pi*sqrt( a^3/mu_S ); % Orbital period [1/s]
tspan = linspace( 0, T, 1000 );
y0_1 = [car_Earth_dep; v_Earth_dep];
[ t_1t, y_1t ] = ode45(@(t_1t, y_1t) ode_2bp(t_1t, y_1t, mu_S), tspan, y0_1, options);

% % 5 Earth departure window
% T = (t12-t11)*24*3600; % duration of the departure window
% tspan = linspace( 0, T, 100 );
% y0_1 = [car_11; v_11];
% [ t_111, y_111 ] = ode45(@(t_111, y_111) ode_2bp(t_111, y_111, mu), tspan, y0_1, options);
% % first point on the departure window, which started at t_11

% 4 Saturn orbit during transfer
tspan = linspace( 0, tspanGA(y)-tspan1(x), 1000 );
y0_2 = [car_Saturn_dep; v_Saturn_dep];
[ t_2, y_2 ] = ode45(@(t_2, y_2) ode_2bp(t_2, y_2, mu_S), tspan, y0_2, options);

% 5 Saturn orbit
a = kepSaturn_dep(1); % Semi-major axis [km]
T = 2*pi*sqrt( a^3/mu_S ); % Orbital period [1/s]
tspan = linspace( 0, T, 1000 );
y0_2 = [car_Saturn_dep; v_Saturn_dep];
[ t_22, y_22 ] = ode45(@(t_22, y_22) ode_2bp(t_22, y_22, mu_S), tspan, y0_2, options);

% 6 Saturn orbit arrival window
T = (t22-t21)*24*3600;
tspan = linspace( 0, T, 1000 );
y0_222 = [car_21; v_21];
options = odeset( 'RelTol', 1e-13, 'AbsTol', 1e-14 );
[ t_222, y_222 ] = ode45(@(t_222, y_222) ode_2bp(t_222, y_222, mu_S), tspan, y0_222, options);


%% PLOT
clc

figure()
hold on
plot3( y_t1(:,1), y_t1(:,2), y_t1(:,3), '-g') % arc of transfer orbit
plot3( y_1t(:,1), y_1t(:,2), y_1t(:,3), '--b') % total Earth orbit
plot3( y_1(:,1), y_1(:,2), y_1(:,3), '-b') % Earth motion during transfer
plot3( y_22(:,1), y_22(:,2), y_22(:,3), '--r') %total Saturn orbit
plot3( y_2(:,1), y_2(:,2), y_2(:,3), '-r') %Saturn motion during transfer
% h1=plot3( y_111(:,1), y_111(:,2), y_111(:,3),'-b', 'LineWidth',20); %departure window
h1.Color(4) = 0.25;
h2=plot3( y_222(:,1), y_222(:,2), y_222(:,3),'-r','LineWidth',20); %arrival window
h2.Color(4) = 0.25;

% drawing planets
scatter3(car_Earth_dep(1), car_Earth_dep(2), car_Earth_dep(3), 100, 'filled', 'b') %optimal point of Earth orbit
scatter3(car_Saturn_fb(1), car_Saturn_fb(2), car_Saturn_fb(3), 120,  'filled', 'r') %optimal point of Saturn orbit
scatter3(y_1(end,1), y_1(end,2), y_1(end,3), 90, 'filled', 'b') % final point of the Earth motion during transfer
scatter3(y_2(end,1), y_2(end,2), y_2(end,3), 120,  'filled', 'g') % final point of the Saturn motion during transfer
scatter3(0, 0, 0, 1000,  'filled', 'y') % Sun

xlabel('X [km]'); ylabel('Y [km]'); zlabel('Z [km]');
title('Transfer problem');
legend('transfer arc', 'Earth orbit','Earth orbit during transfer', 'Saturn orbit', 'Saturn orbit during transfer','Arrival window',...
    '','','','','');
axis equal;
grid on;