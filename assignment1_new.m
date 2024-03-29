clc
clear
close all

% constants
G=astroConstants(1);
mu_S = astroConstants(4);
mu_Saturn = astroConstants(16);
mu_E= astroConstants(13);
R_E= astroConstants(23);
AU= astroConstants(2);
R_Saturn=astroConstants(26);
% Computation of masses, using mu = mass * G
m_Saturn=mu_Saturn/G;
m_Sun=mu_S/G;
m_E=mu_E/G;

% data
dept1=[2032,01,30,0,0,0]; %GIVEN
arrt2=[2062,07,28,0,0,0]; % GIVEN latest arrival on the asteroid
% unused

rp_min=R_E*21;
%% time windows

%synodic period
a_Saturn=9.5826*AU;
T_Saturn = 2*pi*sqrt( a_Saturn^3/mu_S ); % Orbital period [1/s]
a_Earth=AU;
T_Earth=2*pi*sqrt( a_Earth^3/mu_S );
Tsyn_ES=T_Saturn*T_Earth/(abs(T_Saturn-T_Earth));

% Time to go to Saturn from Earth is around 8 years:
%arrt1=[2048,07,28,0,0,0];  %First possible arrival time at the asteroid
%GA_window_1=[2037,01,30,0,0,0]; %assumption for fly-by window
%GA_window_2=[2042,01,30,0,0,0]; %assumption for fly-by window


%[kepNEO_arr,~,~] = ephNEO(date2mjd2000(GA_window_2),86);
[kepNEO_arr,~,~] = ephNEO(0,86);
T_NEO=2*pi*sqrt( (kepNEO_arr(1))^3/mu_S ); % Orbital period [1/s]
Tsyn_SN=T_Saturn*T_NEO/(abs(T_Saturn-T_NEO));
% GA_window_2_mjd2000 = date2mjd2000(GA_window_2)+10*Tsyn_SN/3600/24;
%arrt1=mjd20002date(dept2_mjd2000); %%THIS WAS GIVING US AN ERROR

N_synodic_T=5;

% end of departure window

dept2_mjd2000 = date2mjd2000(dept1)+N_synodic_T*Tsyn_ES/3600/24;
dept2=mjd20002date(dept2_mjd2000);


% Hohmann
a_hohmann=(AU+1.427*10^9)/2;
t_hohmann=pi*sqrt(a_hohmann^3/mu_S);
t_hohmann/3600/24/365.25;

ToF1_min = 6.2; %years
ToF1_max = 10; %years

GA_window_1=dept1 + [ToF1_min, 0, 0 ,0 ,0, 0]; %assumption for fly-by window
GA_window_2=dept2 + [ToF1_max, 0, 0 ,0 ,0, 0]; %assumption for fly-by window
% GA_window_2=dept2 + [ToF1_max, 0, 0 ,0 ,0, 0]; %assumption for fly-by window


% ToF2_min = (1.2*Tsyn_SN)/3600/24; % was 1.2
% ToF2_min_years=ToF2_min/365
% ToF2_max = (3*Tsyn_SN)/3600/24;
% ToF2_max_years=ToF2_max/365

ToF2_min=7; % years to days
ToF2_max = 14;

arrt1 = mjd20002date(date2mjd2000(GA_window_1)+ToF2_min); %First possible arrival time at the asteroid
arrt2 = mjd20002date(date2mjd2000(GA_window_2)+ToF2_max);
%arrt2 = GA_window_2 + [10, 0, 0 ,0 ,0, 0];

% % hypothetic interplanetary parabola with p on the orbit of Earth
% % (increased by a certain factor) and intersecting the orbit of Saturn at
% % theta_parabolic
% p_parabolic=a_Earth*1.3
% theta_parabolic=acos(2*p_parabolic/a_Saturn-1);
% th_deg=theta_parabolic*180/pi
% D=tan(theta_parabolic/2);
% M=1/2*(D+D^3/3);
% n_parabolic=sqrt(mu_S/p_parabolic^3);
% t_parabolic=M/n_parabolic
% t_parabolic/3600/24


%conversion to Modern Julian Date 2000
t_dept1 = date2mjd2000(dept1);
t_dept2 = date2mjd2000(dept2);
t_arrt1 = date2mjd2000(arrt1);
t_arrt2 = date2mjd2000(arrt2);
t_GA_window_1= date2mjd2000(GA_window_1);
t_GA_window_2= date2mjd2000(GA_window_2);

tspan_dept= linspace(t_dept1, t_dept2, 100); %departure window
tspan_arrt = linspace(t_arrt1, t_arrt2, 100); %arrival window
% was 500
tspan_GA=linspace(t_GA_window_1, t_GA_window_2, 100); % fly-by window

p1=3; %Earth
p2=6; %Saturn

%% grid search
ToF1_vect=linspace(ToF1_min*365.25,ToF1_max*365.25,100);
ToF2_vect=linspace(ToF2_min*365.25,ToF2_max*365.25,100);
for i=1:length(tspan_dept)
    i
%     tGA=tspan_dept(i)+ToF1_max;
%     for j=1:length(tspan_GA)       
        for j=1:length(ToF1_vect)
        [dv_1(i,j),V_SC_Saturn_1,V_Saturn, r_Saturn,ToF1(i,j),tpar1(i,j)] = dv_arc1(tspan_dept(i), tspan_dept(i)+ToF1_vect(j), p1, p2, mu_S);        
        % dv_1 manoeuver at Earth
        V_per_min(i,j,:)=V_SC_Saturn_1;       
        for k=1:length(ToF2_vect)
        %for k=1:length(tspan_arrt)             
            [kepNEO_arr,~,~] = ephNEO(tspan_dept(i)+ToF1_vect(j)+ToF2_vect(k),86);
            [dv_2(i,j,k),V_SC_Saturn_2,ToF2(i,j,k),tpar2(i,j,k)] = dv_arc2(tspan_dept(i)+ToF1_vect(j), tspan_dept(i)+ToF1_vect(j)+ToF2_vect(k), r_Saturn, kepNEO_arr, mu_S);
            [rp(i,j,k), Delta_vp(i,j,k),vp_minus(i,j,k),vp_plus(i,j,k),v_inf_minus(i,j,k,:),v_inf_plus(i,j,k,:)] = PGA (V_Saturn, V_SC_Saturn_1',V_SC_Saturn_2', rp_min,mu_Saturn);
            V_per_plus(i,j,k,:) = V_SC_Saturn_2;
            %[rp, Delta_vp] = PGA (V_P,V_minus,V_plus,rp_min,mu_E)
            dv_tot(i,j,k) = dv_1(i,j) + dv_2(i,j,k) + Delta_vp(i,j,k);
        end
    end
end

m1=min(dv_tot);
m2=min(min(dv_tot));
m3=min(min(min(dv_tot)));

% [x,y,z]=find(dv_tot==m3);
[x,y,z] = ind2sub(size(dv_tot),find(dv_tot==m3));

rp_opt=rp(x,y,z);
%% mission results
optimal_dv1=dv_1(x,y)
optimal_dv2=dv_2(x,y,z)
optimal_dvp=Delta_vp(x,y,z)
optimal_departure=mjd20002date(tspan_dept(x))
optimal_Saturn_arrival=mjd20002date(tspan_dept(x)+ToF1_vect(y))
optimal_NEO_arrival=mjd20002date(tspan_dept(x)+ToF1_vect(y)+ToF2_vect(z))

ToF1_vect(y)/365.25
ToF2_vect(z)/365.25


% tpar1(x,y)/3600/24
% tpar2(x,y,z)/3600/24/365

%% Optimal solution - Fly-By time
T=3600*24*2;

v_inf_min=squeeze(v_inf_minus(x,y,z,:));
v_inf_pl=squeeze(v_inf_plus(x,y,z,:));
ecc_minus= 1+(rp(x,y,z)*(norm(v_inf_min)^2)/mu_Saturn);
ecc_plus=1+(rp(x,y,z)*(norm(v_inf_pl)^2)/mu_Saturn);
a_minus=rp(x,y,z)/(1-ecc_minus);
a_plus=rp(x,y,z)/(1-ecc_plus);

% radius of SOI:
r_SOI= a_Saturn*(m_Saturn/m_Sun)^(2/5);
theta_minus= acos((1/ecc_minus)*((a_minus*(1-ecc_minus^2)/r_SOI)-1));
E_minus=2*atanh(sqrt((ecc_minus-1)/(1+ecc_minus)*tan(theta_minus/2)));

theta_plus= acos((1/ecc_plus)*((a_plus*(1-ecc_plus^2)/r_SOI)-1));
E_plus=2*atanh(sqrt((ecc_plus-1)/(1+ecc_plus)*tan(theta_plus/2)));

% Hyperbolic time law:
n_minus=sqrt(mu_Saturn/abs(a_minus)^3);
n_plus=sqrt(mu_Saturn/abs(a_plus)^3);
Delta_t_minus=(ecc_minus*sinh(E_minus)-E_minus)/n_minus;
Delta_t_plus=(ecc_plus*sinh(E_plus)-E_plus)/n_plus;
Delta_t_FlyBy=Delta_t_minus+ Delta_t_plus;
Delta_t_FlyBy_days=Delta_t_FlyBy/3600/24;

%% Porkchop plot Saturn
clc
figure
[X, Y] = meshgrid(tspan_dept, tspan_dept+ToF1_vect);
Z = dv_porkchop(X, Y, p1, p2, @dv_arc1,mu_S);

V=1:2:20;
contour(X./365.25 + 2000 , Y./365.25 + 2000, Z, V,'ShowText','on');
grid on
grid minor
c_porkchop_p=colorbar;
c_porkchop_p.Label.String = ' Δv [s]';

xlabel('Departure time [years]');
ylabel('Arrival time [years]');
title('Transfer to Saturn Porkchop plot');
hold on;
scatter(tspan_dept(x)/365.25+2000,(tspan_dept(x)+ToF1_vect(y))/365.25+2000,25,'red','filled');

% figure
% grid off
% surface(X./365.25 + 2000,Y./365.25 + 2000,Z)
% colorbar
% xlabel('Departure time [years]');
% ylabel('Arrival time [years]');
% hold on
% scatter(tspan_dept(x)/365.25+2000,tspan_GA(y)/365.25+2000,300,'red','filled');

%% Porkchop plot asteroid
figure()

[X, Y] = meshgrid(tspan_dept+ToF1_vect, tspan_dept+ToF1_vect+ToF2_vect);
Z = dv_porkchop(X, Y, p2,86, @dv_arcNEO,mu_S);

V=2:30;
contour(X./365.25 +2000 , Y./365.25 + 2000, Z, V,'ShowText','on');
grid on
grid minor
colorbar
xlabel('Departure time [years]');
ylabel('Arrival time [years]');
title('Transfer to Asteroid Porkchop plot');
hold on;
scatter((tspan_dept(x)+ToF1_vect(y))/365.25+2000,(tspan_dept(x)+ToF1_vect(y)+ToF2_vect(z))/365.25+2000,25,'red','filled');

% %% cost plot ????
% clc
% NEO=86;
% % tspan_1=linspace(t_dept1,t_GA_window_2, 1000);
% % tspan_2=linspace(t_GA_window_2,t_arrt2, 1000);
% 
% % [X, Y] = meshgrid(tspan_1, tspan_2);
% 
% % [X, Y] = meshgrid(tspan_dept, tspan_arrt);
%         % Saturn    
%         [X1, Y1_GA] = meshgrid(tspan_dept, tspan_GA);
%         % Z = dv_porkchop(X, Y, p1, p2, @dv_arc1,mu_S);
% 
%         % Asteroid
%         [X2_GA, Y2] = meshgrid(tspan_GA, tspan_arrt);
%         % Z = dv_porkchop(X, Y, p2,86, @dv_arcNEO,mu_S);
% [Z,X,Y] = dv_totalcost(X1, Y1_GA,X2_GA,Y2, p1, p2,NEO, mu_S);
% % function Z = dv_totalcost(X1, Y1_GA,X2_GA,Y2, p1, p2,NEO,mu)
% % X=[X1,X2_GA]
% % Y=[Y1_GA;Y2]
% V=10:2:40;
% contour(X./365.25 +2000 , Y./365.25 + 2000, Z, V,'ShowText','on');
% colorbar
% xlabel('Departure time [years]');
% ylabel('Arrival time [years]');
% title('Cost plot');
% hold on;
% % scatter(tspan_dept(x)/365.25+2000,tspan_arrt(z)/365.25+2000,20,'red','filled');

%% Plot the planetocentric hyperbolic arcs:
figure()
earth_sphere
hold on
grid on

y0 =[rp(x,y,z)*[1;0;0]; vp_minus(x,y,z)*[0;1;0]];

tspan = linspace( 0, -2*T,1000);
options = odeset( 'RelTol', 1e-14, 'AbsTol', 1e-14 );
[t, Y_FlyBy_before ] = ode113( @(t,y) ode_2bp(t,y,mu_Saturn), tspan, y0, options);
plot3( Y_FlyBy_before(:,1), Y_FlyBy_before(:,2), Y_FlyBy_before(:,3), '-b','LineWidth',2);
hold on

y0=[rp(x,y,z)*[1;0;0];vp_plus(x,y,z)*[0;1;0]];
tspan = linspace( 0, 2*T,1000);
options = odeset( 'RelTol', 1e-14, 'AbsTol', 1e-14 );
[t, Y_FlyBy_after ] = ode113( @(t,y) ode_2bp(t,y,mu_Saturn), tspan, y0, options);
plot3( Y_FlyBy_after(:,1), Y_FlyBy_after(:,2), Y_FlyBy_after(:,3), '-r','LineWidth',2);
xlabel('x [ R_{Saturn} ]');
ylabel('y [ R_{Saturn} ]');
zlabel('z [ R_{Saturn} ]');
hold on

%Entry Asymptote
x_asymp1= -18 *R_Saturn  : rp(x,y,z)-a_minus;
y_asymp1= -tan(acos(-1/ecc_minus))*(x_asymp1+a_minus-rp(x,y,z));
plot(x_asymp1,y_asymp1,'b-',LineWidth=1);

%Exit Asymptote
x_asymp2=-18*R_Saturn:rp(x,y,z)-a_plus;
y_asymp2= tan(acos(-1/ecc_plus))*(x_asymp2+a_plus-rp(x,y,z));
plot(x_asymp2,y_asymp2,'r-',LineWidth=1);

x_3=-60*R_Saturn:30*R_Saturn;
y_3=0*x_3;
plot(x_3,y_3,'--k',LineWidth=1);
%  axis([-20*R_Saturn 20*R_Saturn -20*R_Saturn 20*R_Saturn]);
scatter3(0,0,0,80,'red','filled');
legend('','Hyperbola entry leg', 'Hyperbola exit leg','Entry Asymptote', 'Exit Asymptote','','Saturn');
%% heliocentric plot

figure()
earth_sphere('AU')
hold on
%Plotting of 1st transfer arc
[kep_1,~] = uplanet(tspan_dept(x)+ToF1_vect(y),p2);
[r0_Saturn,v0_Saturn] = par2car(kep_1(1),kep_1(2),kep_1(3),kep_1(4),kep_1(5),kep_1(6),mu_S);
V_minus = squeeze(V_per_min(x,y,:));
% V_minus = V_per_min;
y0_SC_saturn=[r0_Saturn;V_minus];
% T=3600*24*365*5;
% % T=2*pi*sqrt( a^3/mu_E ); % Orbital period [1/s]

% tspan = linspace( 0, -(tspan_GA(y)-tspan_dept(x))*24*3600,1000);
tspan = linspace( 0, -(ToF1_vect(y))*24*3600,1000);

     % Set options for the ODE solver
options = odeset( 'RelTol', 1e-14, 'AbsTol', 1e-14 );

%Heliocentric transfer arc before FlyBy:
[ t, Y_helio_before ] = ode113( @(t,y) ode_2bp(t,y,mu_S), tspan, y0_SC_saturn, options );
plot3( Y_helio_before(:,1)/AU, Y_helio_before(:,2)/AU, Y_helio_before(:,3)/AU, '-','LineWidth',2);
hold on
grid on
%Sun sphere
scatter3(0, 0 ,0 ,100, 'yellow', 'filled');
%Saturn sphere
scatter3(r0_Saturn(1)/AU,r0_Saturn(2)/AU,r0_Saturn(3)/AU,20,'red','filled');

%Plotting of the 2nd transfer arc:
V_plus = squeeze(V_per_plus(x,y,z,:));
y0_SC_saturn_plus=[r0_Saturn;V_plus];
% T=2*pi*sqrt( a^3/mu_E ); % Orbital period [1/s]
tspan = linspace( 0, (ToF2_vect(z))*24*3600,1000);
     % Set options for the ODE solver
options = odeset( 'RelTol', 1e-14, 'AbsTol', 1e-14 );
[ t, Y_helio_after ] = ode113( @(t,y) ode_2bp(t,y,mu_S), tspan, y0_SC_saturn_plus, options );
plot3( Y_helio_after(:,1)/AU, Y_helio_after(:,2)/AU, Y_helio_after(:,3)/AU, 'r-','LineWidth',2);
legend('','Transfer arc 1','Sun','Saturn','Transfer arc 2','Saturn orbit','Earth orbit','Earth','Asteroid orbit','Asteroid');

%Plotting of objects orbits

%Saturn orbit:
T=2*pi*sqrt(kep_1(1)^3/mu_S);
tspan = linspace( 0, T,1000);
y0_Saturn = [r0_Saturn,v0_Saturn];
     % Set options for the ODE solver
options = odeset( 'RelTol', 1e-14, 'AbsTol', 1e-14 );
[ t, Y_saturn ] = ode113( @(t,y) ode_2bp(t,y,mu_S), tspan, y0_Saturn, options );
plot3( Y_saturn(:,1)/AU, Y_saturn(:,2)/AU, Y_saturn(:,3)/AU, 'r--','LineWidth',2);

%Earth orbit
[kep_2,~] = uplanet(tspan_dept(x),p1); %Time doesn't matter, because we do it for one period
[r0_Earth,v0_Earth] = par2car(kep_2(1),kep_2(2),kep_2(3),kep_2(4),kep_2(5),kep_2(6),mu_S);
y0_Earth = [r0_Earth;v0_Earth];
T=2*pi*sqrt(kep_2(1)^3/mu_S);
tspan = linspace( 0, T,1000);
     % Set options for the ODE solver
options = odeset( 'RelTol', 1e-14, 'AbsTol', 1e-14 );
[ t, Y_earth ] = ode113( @(t,y) ode_2bp(t,y,mu_S), tspan, y0_Earth, options );
plot3( Y_earth(:,1)/AU, Y_earth(:,2)/AU, Y_earth(:,3)/AU, 'b--','LineWidth',2);
scatter3(r0_Earth(1)/AU,r0_Earth(2)/AU,r0_Earth(3)/AU,30,'blue','filled');

%Asteroid orbit
[kep_3,~] = ephNEO(tspan_dept(x),86); %Time doesn't matter, because we do it for one period
[r0_NEO,v0_NEO] = par2car(kep_3(1),kep_3(2),kep_3(3),kep_3(4),kep_3(5),kep_3(6),mu_S);
y0_NEO = [r0_NEO;v0_NEO];
T=(ToF1_vect(y)+ToF2_vect(z))*24*3600; %Converted to seconds
tspan = linspace( 0, T,1000);
     % Set options for the ODE solver
options = odeset( 'RelTol', 1e-14, 'AbsTol', 1e-14 );
[ t, Y_asteroid] = ode113( @(t,y) ode_2bp(t,y,mu_S), tspan, y0_NEO, options );
plot3( Y_asteroid(:,1)/AU, Y_asteroid(:,2)/AU, Y_asteroid(:,3)/AU, '--g','LineWidth',2);
%Getting new ephemerides for the scatter
[kep_3,~] = ephNEO(tspan_dept(x)+ToF1_vect(y)+ToF2_vect(z),86); %Time doesn't matter, because we do it for one period
[r0,v0] = par2car(kep_3(1),kep_3(2),kep_3(3),kep_3(4),kep_3(5),kep_3(6),mu_S);

scatter3(r0_NEO(1)/AU,r0_NEO(2)/AU,r0_NEO(3)/AU,20,'green','filled');

scatter3(Y_earth(end,1)/AU, Y_earth(end,2)/AU, Y_earth(end,3)/AU, 30, 'filled', 'b') % final point of the Earth motion during transfer
scatter3(r0_Saturn(1)/AU,r0_Saturn(2)/AU,r0_Saturn(3)/AU,30,'red','filled'); %Saturn position at departure
scatter3(r0_NEO(1)/AU,r0_NEO(2)/AU,r0_NEO(3)/AU,20,'green','filled'); %Saturn position at departure

legend('','Transfer arc 1','Sun','Saturn','Transfer arc 2','Saturn orbit','Earth orbit','Earth','Asteroid orbit','Asteroid','','','');

%% without arcs

figure()
earth_sphere('AU')
hold on
grid on
scatter3(0, 0 ,0 ,30, 'yellow', 'filled');

%Plotting of 1st transfer arc
[kep_1,~] = uplanet(tspan_dept(x)+ToF1_vect(y),p2);
[r0_Saturn,v0_Saturn] = par2car(kep_1(1),kep_1(2),kep_1(3),kep_1(4),kep_1(5),kep_1(6),mu_S);
V_minus = squeeze(V_per_min(x,y,:));
% V_minus = V_per_min;
y0_SC_saturn=[r0_Saturn;V_minus];
% T=3600*24*365*5;
% % T=2*pi*sqrt( a^3/mu_E ); % Orbital period [1/s]

% tspan = linspace( 0, -(tspan_GA(y)-tspan_dept(x))*24*3600,1000);
tspan = linspace( 0, -(ToF1_vect(y))*24*3600,1000);

     % Set options for the ODE solver
options = odeset( 'RelTol', 1e-14, 'AbsTol', 1e-14 );
%Plotting of objects orbits

%Saturn orbit:
T=2*pi*sqrt(kep_1(1)^3/mu_S);
tspan = linspace( 0, T,1000);
y0_Saturn = [r0_Saturn,v0_Saturn];
     % Set options for the ODE solver
options = odeset( 'RelTol', 1e-14, 'AbsTol', 1e-14 );
[ t, Y_saturn ] = ode113( @(t,y) ode_2bp(t,y,mu_S), tspan, y0_Saturn, options );
plot3( Y_saturn(:,1)/AU, Y_saturn(:,2)/AU, Y_saturn(:,3)/AU, 'r','LineWidth',2);

%Earth orbit
[kep_2,~] = uplanet(tspan_dept(x),p1); %Time doesn't matter, because we do it for one period
[r0_Earth,v0_Earth] = par2car(kep_2(1),kep_2(2),kep_2(3),kep_2(4),kep_2(5),kep_2(6),mu_S);
y0_Earth = [r0_Earth;v0_Earth];
T=2*pi*sqrt(kep_2(1)^3/mu_S);
tspan = linspace( 0, T,1000);
     % Set options for the ODE solver
options = odeset( 'RelTol', 1e-14, 'AbsTol', 1e-14 );
[ t, Y_earth ] = ode113( @(t,y) ode_2bp(t,y,mu_S), tspan, y0_Earth, options );
plot3( Y_earth(:,1)/AU, Y_earth(:,2)/AU, Y_earth(:,3)/AU, 'b','LineWidth',2);

%Asteroid orbit
[kep_3,~] = ephNEO(tspan_dept(x),86); %Time doesn't matter, because we do it for one period
[r0_NEO,v0_NEO] = par2car(kep_3(1),kep_3(2),kep_3(3),kep_3(4),kep_3(5),kep_3(6),mu_S);
y0_NEO = [r0_NEO;v0_NEO];
T=(ToF1_vect(y)+ToF2_vect(z))*24*3600; %Converted to seconds
tspan = linspace( 0, T,1000);
     % Set options for the ODE solver
options = odeset( 'RelTol', 1e-14, 'AbsTol', 1e-14 );
[ t, Y_asteroid] = ode113( @(t,y) ode_2bp(t,y,mu_S), tspan, y0_NEO, options );
plot3( Y_asteroid(:,1)/AU, Y_asteroid(:,2)/AU, Y_asteroid(:,3)/AU, 'g','LineWidth',2);
%Getting new ephemerides for the scatter
[kep_3,~] = ephNEO(tspan_dept(x)+ToF1_vect(y)+ToF2_vect(z),86); %Time doesn't matter, because we do it for one period
[r0,v0] = par2car(kep_3(1),kep_3(2),kep_3(3),kep_3(4),kep_3(5),kep_3(6),mu_S);

legend('','Sun','Saturn orbit','Earth orbit','Asteroid orbit');
