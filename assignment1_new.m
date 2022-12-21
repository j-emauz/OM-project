clear
clc
close all

% constants
mu_S = astroConstants(4);
mu_Saturn = astroConstants(16);
mu_E= astroConstants(13);
R_E= astroConstants(23);
AU= astroConstants(2);

% data
dept1=[2028,01,30,0,0,0]; %GIVEN
arrt2=[2062,07,28,0,0,0]; % GIVEN latest arrival on the asteroid
rp_min=R_E*21 * 1.2;

% synodic period
a_Saturn=9.5826*AU;
T_Saturn = 2*pi*sqrt( a_Saturn^3/mu_S ); % Orbital period [1/s]
a_Earth=AU;
T_Earth=2*pi*sqrt( a_Earth^3/mu_S );
Tsyn_ES=T_Saturn*T_Earth/(abs(T_Saturn-T_Earth));

N_synodic_T=5;
% end of departure window
dept2_mjd2000 = date2mjd2000(dept1)+N_synodic_T*Tsyn_ES/3600/24;
dept2=mjd20002date(dept2_mjd2000);

% theta_parabolic=2/3*pi;
% p_parabolic=100;
% D=tan(theta_parabolic/2);
% M=1/2*(D+D^3/3);
% n_parabolic=sqrt(mu_S/p_parabolic^3);
% t_parabolic=M/n_parabolic;
% dept2_mjd2000 = date2mjd2000(dept1)

% Time to go to Saturn from Earth is around 8 years:
arrt1=[2048,07,28,0,0,0];  %First possible arrival time at the asteroid
GA_window_1=[2037,01,30,0,0,0]; %assumption for fly-by window
GA_window_2=[2042,01,30,0,0,0]; %assumption for fly-by window

[kepNEO_arr,~,~] = ephNEO(date2mjd2000(GA_window_2),86);
T_NEO=2*pi*sqrt( (kepNEO_arr(1))^3/mu_S ); % Orbital period [1/s]
Tsyn_SN=T_Saturn*T_NEO/(abs(T_Saturn-T_NEO));
GA_window_2_mjd2000 = date2mjd2000(GA_window_2)+10*Tsyn_SN/3600/24;
%arrt1=mjd20002date(dept2_mjd2000); %%THIS WAS GIVING US AN ERROR

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


for i=1:length(tspan_dept)
    i
    for j=1:length(tspan_GA)       
        [dv_1(i,j),V_SC_Saturn_1,V_Saturn, r_Saturn, t1(i,j),ToF1] = dv_arc1(tspan_dept(i), tspan_GA(j), p1, p2, mu_S);        
        % dv_1 manoeuver at Earth
        V_per_min(i,j,:)=V_SC_Saturn_1;
        for k=1:length(tspan_arrt)             
            [kepNEO_arr,~,~] = ephNEO(tspan_arrt(k),86);
            [dv_2(i,j,k),V_SC_Saturn_2, t2(i,j,k),ToF2] = dv_arc2(tspan_GA(j), tspan_arrt(k), r_Saturn, kepNEO_arr, mu_S);
            [rp, Delta_vp(i,j,k)] = PGA (V_Saturn, V_SC_Saturn_1',V_SC_Saturn_2', rp_min,mu_Saturn);
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

%% Porkchop plot Saturn
clc

[X, Y] = meshgrid(tspan_dept, tspan_GA);
Z = dv_grid_arc1(X, Y, p1, p2, mu_S);

V=5:30;
contour(X, Y, Z, V,'ShowText','on');
colorbar
% surface(X,Y,Z)

%% Porkchop plot asteroid

[X, Y] = meshgrid(tspan_GA, tspan_arrt);
Z = dv_calcgrid(X, Y, p1, p2, mu_S);

V=5:30;
contour(X, Y, Z, V,'ShowText','on');
colorbar
%% Plot the planetocentric hyperbolic arcs:
earth_sphere
hold on
grid on
T=3600*2;
%y0 =
tspan = linspace( 0, -T,1000);
options = odeset( 'RelTol', 1e-14, 'AbsTol', 1e-14 );
[t, Y_planet_before ] = ode113( @(t,y) ode_2bp(t,y,mu_E), tspan, y0, options);
plot3( Y_planet_before(:,1), Y_planet_before(:,2), Y_planet_before(:,3), '-','LineWidth',2);
hold on

T=3600*2;
y0=[rp*[1;0;0]; vp_plus*[0;1;0]];
tspan = linspace( 0, T,1000);
options = odeset( 'RelTol', 1e-14, 'AbsTol', 1e-14 );
[t, Y_planet_before ] = ode113( @(t,y) ode_2bp(t,y,mu_E), tspan, y0, options);
plot3( Y_planet_before(:,1), Y_planet_before(:,2), Y_planet_before(:,3), '-','LineWidth',2);
xlabel('x [R_E]');
ylabel('y [R_E]');
zlabel('z [R_E]');
hold on

x= -4 *R_E  : rp-a_minus;
y= -tan(acos(-1/ecc_minus(rp)))*(x+a_minus-rp);
plot(x,y,'b-',LineWidth=2);

hold on
x=-4*R_E:rp-a_plus;
y= tan(acos(-1/ecc_plus(rp)))*(x+a_plus-rp);
plot(x,y,'r-',LineWidth=2);

x=-6*R_E:6*R_E;
y=0*x;
plot(x,y,'--k',LineWidth=2);
axis([-5*R_E 5*R_E -10*R_E 10*R_E]);
%% heliocentric leg

figure()
earth_sphere('AU')
hold on
%Plotting of 1st transfer arc
[kep_1,~] = uplanet(tspan_GA(y),p2);
[r0,v0] = par2car(kep_1(1),kep_1(2),kep_1(3),kep_1(4),kep_1(5),kep_1(6),mu_S);
%V_minus = V_per_min(x,y,:);
V_minus = V_per_min;
y0=[r0;V_minus];
% T=3600*24*365*5;
% % T=2*pi*sqrt( a^3/mu_E ); % Orbital period [1/s]
tspan = linspace( 0, -(tspan_GA(y)-tspan_dept(x))*24*3600,1000);
     % Set options for the ODE solver
options = odeset( 'RelTol', 1e-14, 'AbsTol', 1e-14 );
[ t, Y_helio_before ] = ode113( @(t,y) ode_2bp(t,y,mu_S), tspan, y0, options );
plot3( Y_helio_before(:,1)/AU, Y_helio_before(:,2)/AU, Y_helio_before(:,3)/AU, '-','LineWidth',2);
hold on
grid on
%Sun sphere
scatter3(0, 0 ,0 ,100, 'yellow', 'filled');
%Saturn sphere
scatter3(r0(1)/AU,r0(2)/AU,r0(3)/AU,20,'red','filled');

%Plotting of the 2nd transfer arc
V_plus = squeeze(V_per_plus(x,y,z,:));
y0=[r0;V_plus];
% T=2*pi*sqrt( a^3/mu_E ); % Orbital period [1/s]
tspan = linspace( 0, (tspan_arrt(z)-tspan_GA(y))*24*3600,1000);
     % Set options for the ODE solver
options = odeset( 'RelTol', 1e-14, 'AbsTol', 1e-14 );
[ t, Y_helio_before ] = ode113( @(t,y) ode_2bp(t,y,mu_S), tspan, y0, options );
plot3( Y_helio_before(:,1)/AU, Y_helio_before(:,2)/AU, Y_helio_before(:,3)/AU, 'r-','LineWidth',2);
legend('','Transfer arc 1','Sun','Saturn','Transfer arc 2','Saturn orbit','Earth orbit','Earth','Asteroid orbit','Asteroid');

%Plotting of objects orbits

%Saturn orbit
T=2*pi*sqrt(kep_1(1)^3/mu_S);
tspan = linspace( 0, T,1000);
y0 = [r0,v0];
     % Set options for the ODE solver
options = odeset( 'RelTol', 1e-14, 'AbsTol', 1e-14 );
[ t, Y_saturn ] = ode113( @(t,y) ode_2bp(t,y,mu_S), tspan, y0, options );
plot3( Y_saturn(:,1)/AU, Y_saturn(:,2)/AU, Y_saturn(:,3)/AU, 'r--','LineWidth',2);

%Earth orbit
[kep_2,~] = uplanet(tspan_dept(x),p1); %Time doesn't matter, because we do it for one period
[r0,v0] = par2car(kep_2(1),kep_2(2),kep_2(3),kep_2(4),kep_2(5),kep_2(6),mu_S);
y0 = [r0;v0];
T=2*pi*sqrt(kep_2(1)^3/mu_S);
tspan = linspace( 0, T,1000);
     % Set options for the ODE solver
options = odeset( 'RelTol', 1e-14, 'AbsTol', 1e-14 );
[ t, Y_earth ] = ode113( @(t,y) ode_2bp(t,y,mu_S), tspan, y0, options );
plot3( Y_earth(:,1)/AU, Y_earth(:,2)/AU, Y_earth(:,3)/AU, 'b--','LineWidth',2);
scatter3(r0(1)/AU,r0(2)/AU,r0(3)/AU,30,'blue','filled');

%Asteroid orbit
[kep_3,~] = ephNEO(tspan_dept(x),86); %Time doesn't matter, because we do it for one period
[r0,v0] = par2car(kep_3(1),kep_3(2),kep_3(3),kep_3(4),kep_3(5),kep_3(6),mu_S);
y0 = [r0;v0];
T=(tspan_arrt(z) - tspan_dept(x))*24*3600; %Converted to seconds
tspan = linspace( 0, T,1000);
     % Set options for the ODE solver
options = odeset( 'RelTol', 1e-14, 'AbsTol', 1e-14 );
[ t, Y_asteroid] = ode113( @(t,y) ode_2bp(t,y,mu_S), tspan, y0, options );
plot3( Y_asteroid(:,1)/AU, Y_asteroid(:,2)/AU, Y_asteroid(:,3)/AU, 'g--','LineWidth',2);
%Getting new ephemerides for the scatter
[kep_3,~] = ephNEO(tspan_arrt(z),86); %Time doesn't matter, because we do it for one period
[r0,v0] = par2car(kep_3(1),kep_3(2),kep_3(3),kep_3(4),kep_3(5),kep_3(6),mu_S);
scatter3(r0(1)/AU,r0(2)/AU,r0(3)/AU,20,'green','filled');
legend('','Transfer arc 1','Sun','Saturn','Transfer arc 2','Saturn orbit','Earth orbit','Earth','Asteroid orbit','Asteroid');