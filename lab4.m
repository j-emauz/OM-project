% 1.a
clc
clear
close all

% constants
mu_S = astroConstants(4);
mu_E= astroConstants(13);
AU = astroConstants(2);
R_E=astroConstants(23);

% data
v_inf_minus = [15.1; 0; 0]; % km/s
v_inf = norm(v_inf_minus);
Delta = 9200; %km 
r_E= [1; 0; 0]*AU; % AU

up = [0; 0; -1]; %under planet
fp = [0; 1; 0]; %front planet
bp = [0; -1; 0]; %behind planet


uu = cross(Delta*up, v_inf_minus)/norm(cross(Delta*up, v_inf_minus));
uf = cross(Delta*fp, v_inf_minus)/norm(cross(Delta*fp, v_inf_minus));
ub = cross(Delta*bp, v_inf_minus)/norm(cross(Delta*bp, v_inf_minus));

a = -mu_E/(v_inf^2);

delta = 2*atan2(-a, Delta);
e = 1/sin(delta/2);
rp = a*(1-e);

v_inf_plus_u = v_inf_minus*cos(delta) + cross(uu, v_inf_minus)*sin(delta) + uu*(dot(uu, v_inf_minus))*(1-cos(delta));
v_inf_plus_f = v_inf_minus*cos(delta) + cross(uf, v_inf_minus)*sin(delta) + uf*(dot(uf, v_inf_minus))*(1-cos(delta));
v_inf_plus_b = v_inf_minus*cos(delta) + cross(ub, v_inf_minus)*sin(delta) + ub*(dot(ub, v_inf_minus))*(1-cos(delta));

delta_v_u = v_inf_plus_u - v_inf_minus;
delta_v_f = v_inf_plus_f - v_inf_minus;
delta_v_b = v_inf_plus_b - v_inf_minus;

V_P = [0; 1; 0]*sqrt(mu_S/norm(r_E));
V_plus_b = V_P + v_inf_plus_b;
V_plus_f = V_P + v_inf_plus_f;
V_plus_u = V_P + v_inf_plus_u;
V_minus = V_P + v_inf_minus;


%% integrating

% time of integration
T=3600*1.5;

%under
r0=[-5*R_E;0; -Delta];
y0=[r0;v_inf_minus];
% T=2*pi*sqrt( a^3/mu_E ); % Orbital period [1/s]
tspan = linspace( 0, T,1000);
     % Set options for the ODE solver
options = odeset( 'RelTol', 1e-14, 'AbsTol', 1e-14 );
[ Tu, Yu ] = ode113( @(t,y) ode_2bp(t,y,mu_E), tspan, y0, options );

[m  pos_u]=min(vecnorm(Yu(:,1:3),2,2));
rp_u=Yu(pos_u,1:3);

% behind
r0=[-5*R_E; -Delta; 0 ];
y0=[r0; v_inf_minus];
% T=2*pi*sqrt( a^3/mu_E ); % Orbital period [1/s]
tspan = linspace( 0, T,1000);
     % Set options for the ODE solver
options = odeset( 'RelTol', 1e-14, 'AbsTol', 1e-14 );
[ Tb, Yb ] = ode113( @(t,y) ode_2bp(t,y,mu_E), tspan, y0, options );

[m  pos_b]=min(vecnorm(Yb(:,1:3),2,2));
rp_b=Yb(pos_b,1:3);

% in front
r0=[-5*R_E; Delta; 0 ];
y0=[r0;v_inf_minus];
% T=2*pi*sqrt( a^3/mu_E ); % Orbital period [1/s]
tspan = linspace( 0, T,1000);
     % Set options for the ODE solver
options = odeset( 'RelTol', 1e-14, 'AbsTol', 1e-14 );
[ Tf, Yf ] = ode113( @(t,y) ode_2bp(t,y,mu_E), tspan, y0, options );

[m  pos_f]=min(vecnorm(Yf(:,1:3),2,2));
rp_f=Yf(pos_f,1:3);
%% 3d plotting
figure(1)
earth_sphere;
grid on
hold on;
plot3( Yu(:,1), Yu(:,2), Yu(:,3), 'r-','LineWidth',2);

plot3( Yb(:,1), Yb(:,2), Yb(:,3), 'g-', 'LineWidth',2);

plot3( Yf(:,1), Yf(:,2), Yf(:,3), 'b-','LineWidth',2 );
scatter3(rp_u(1), rp_u(2), rp_u(3));
scatter3(rp_b(1), rp_b(2), rp_b(3));
scatter3(rp_f(1), rp_f(2), rp_f(3));
legend('','Flyby under the planet', 'Flyby behind the planet','Flyby in front the planet')

%% heliocentric leg
% time of integration
T=3600*24*365/4;

%before flyby
r0=r_E;
y0=[r0;V_minus];
% T=2*pi*sqrt( a^3/mu_E ); % Orbital period [1/s]
tspan = linspace( 0, -T,1000);
     % Set options for the ODE solver
options = odeset( 'RelTol', 1e-14, 'AbsTol', 1e-14 );
[ t, Y_helio_before ] = ode113( @(t,y) ode_2bp(t,y,mu_S), tspan, y0, options );

% after flyby u
r0=r_E;
y0=[r0;V_plus_u];
% T=2*pi*sqrt( a^3/mu_E ); % Orbital period [1/s]
tspan = linspace( 0, T,1000);
     % Set options for the ODE solver
options = odeset( 'RelTol', 1e-14, 'AbsTol', 1e-14 );
[ t, Y_helio_after ] = ode113( @(t,y) ode_2bp(t,y,mu_S), tspan, y0, options );

% plotting
figure (2)
grid on
axis([-1.5 1.5 -1.5 1.5 -1 1]);
plot3( Y_helio_before(:,1)/AU, Y_helio_before(:,2)/AU, Y_helio_before(:,3)/AU, 'b-','LineWidth',2);
axis([-1.5 1.5 -1.5 1.5 -1 1]);
hold on
plot3( Y_helio_after(:,1)/AU, Y_helio_after(:,2)/AU, Y_helio_after(:,3)/AU, 'r-','LineWidth',2 );
axis([-1.5 1.5 -1.5 1.5 -1 1]);
hold on
scatter3(0, 0 ,0 ,'yellow', 'filled');
axis([-1.5 1.5 -1.5 1.5 -1 1]);


% after flyby b
r0=r_E;
y0=[r0;V_plus_b];
% T=2*pi*sqrt( a^3/mu_E ); % Orbital period [1/s]
tspan = linspace( 0, T,1000);
     % Set options for the ODE solver
options = odeset( 'RelTol', 1e-14, 'AbsTol', 1e-14 );
[ t, Y_helio_after ] = ode113( @(t,y) ode_2bp(t,y,mu_S), tspan, y0, options );

% plotting
figure (3)
grid on
plot3( Y_helio_before(:,1)/AU, Y_helio_before(:,2)/AU, Y_helio_before(:,3)/AU, 'b-','LineWidth',2);
hold on
plot3( Y_helio_after(:,1)/AU, Y_helio_after(:,2)/AU, Y_helio_after(:,3)/AU, 'r-','LineWidth',2 );
hold on
scatter3(0, 0 ,0 ,'yellow', 'filled');

% after flyby f
r0=r_E;
y0=[r0;V_plus_f];
% T=2*pi*sqrt( a^3/mu_E ); % Orbital period [1/s]
tspan = linspace( 0, T,1000);
     % Set options for the ODE solver
options = odeset( 'RelTol', 1e-14, 'AbsTol', 1e-14 );
[ t, Y_helio_after ] = ode113( @(t,y) ode_2bp(t,y,mu_S), tspan, y0, options );

% plotting
figure (4)
grid on
plot3( Y_helio_before(:,1)/AU, Y_helio_before(:,2)/AU, Y_helio_before(:,3)/AU, 'b-','LineWidth',2);
hold on
plot3( Y_helio_after(:,1)/AU, Y_helio_after(:,2)/AU, Y_helio_after(:,3)/AU, 'r-','LineWidth',2 );
hold on
scatter3(0, 0 ,0 ,'yellow', 'filled');

pause(2)

%% % 1.b
clc
clear
close all

% constants
mu_S = astroConstants(4);
mu_E= astroConstants(13);
AU = astroConstants(2);
R_E=astroConstants(23);

% data
v_inf_minus = [15.1; 0; 0]; % km/s
v_inf = norm(v_inf_minus);
r_E= [1; 0; 0]*AU; % AU

fp = [0; 1; 0]; %front planet

a = -mu_E/(v_inf^2);
V_P = [0; 1; 0]*sqrt(mu_S/norm(r_E));

Delta=0:200:40*R_E;
Delta_vect = Delta/R_E;
for i=1:length(Delta)
%Delta(i)=Delta_vect(i)*R_E; 
uf = cross(Delta(i)*fp, v_inf_minus)/norm(cross(Delta(i)*fp, v_inf_minus));
delta(i)= 2*atan2(-a, Delta(i));
e(i)= 1/sin(delta(i)/2);
rp(i) = a*(1-e(i));

v_inf_plus(:,i) = v_inf_minus*cos(delta(i)) + cross(uf, v_inf_minus)*sin(delta(i)) + uf*(dot(uf, v_inf_minus))*(1-cos(delta(i)));

delta_v(:,i) = v_inf_plus(:,i) - v_inf_minus;

V_plus(:,i) = V_P + v_inf_plus(:,i);

V_minus = V_P + v_inf_minus;

end

figure(1)
plot(Delta_vect,rp/R_E)
xlabel('Impact parameter \delta (over R_E)')
ylabel('Flyby minimum altitude (over R_E)')
figure(2)
plot(Delta_vect,delta*180/pi)
xlabel('Impact parameter Δ (over R_E)')
ylabel('Turning angle [deg]')


%% integrating
clc
% time of integration
T=3600*1.5;
options = odeset( 'RelTol', 1e-14, 'AbsTol', 1e-14 );

D=[9200,10200,11200,12200,13200];

figure()
earth_sphere;
legend_entries={};
for i=1:length(D)
% in front
r0=[-5*R_E; D(i); 0 ];
y0=[r0;v_inf_minus];
% T=2*pi*sqrt( a^3/mu_E ); % Orbital period [1/s]
tspan = linspace( 0, T,1000);
     % Set options for the ODE solver
[ Tf, Yf ] = ode113( @(t,y) ode_2bp(t,y,mu_E), tspan, y0, options );

[m  pos_f]=min(vecnorm(Yf(:,1:3),2,2));
rp_f=Yf(pos_f,1:3);

% 3d plotting

grid on
hold on;
plot3( Yf(:,1), Yf(:,2), Yf(:,3), '-','LineWidth',2 );
scatter3(rp_f(1), rp_f(2), rp_f(3));

legend_entries{end+1}=[''];
legend_entries{end+1}=['Δ_',num2str(i)];
end

legend(legend_entries);


%% heliocentric leg
clc
% time of integration
T=3600*24*365/2;
options = odeset( 'RelTol', 1e-14, 'AbsTol', 1e-14 );

D=[9200,10200,11200,12200,13200];

legend_entries={};
figure()

%before flyby
r0=r_E;
y0=[r0;V_minus];
% T=2*pi*sqrt( a^3/mu_E ); % Orbital period [1/s]
tspan = linspace( 0, -T,1000);
     % Set options for the ODE solver
options = odeset( 'RelTol', 1e-14, 'AbsTol', 1e-14 );
[ t, Y_helio_before ] = ode113( @(t,y) ode_2bp(t,y,mu_S), tspan, y0, options );
plot3( Y_helio_before(:,1)/AU, Y_helio_before(:,2)/AU, Y_helio_before(:,3)/AU, '-','LineWidth',2);
hold on
grid on
scatter3(0, 0 ,0 ,100, 'yellow', 'filled');
hold on
legend_entries{1}=['Orbit before flyby'];
legend_entries{2}=[''];

for i=1:length(D)

% after flyby
r0=r_E;
y0=[r0;V_plus(:,find(Delta==D(i)))];
% T=2*pi*sqrt( a^3/mu_E ); % Orbital period [1/s]
tspan = linspace( 0, T,1000);
     % Set options for the ODE solver
options = odeset( 'RelTol', 1e-14, 'AbsTol', 1e-14 );
[ t, Y_helio_after ] = ode113( @(t,y) ode_2bp(t,y,mu_S), tspan, y0, options );

plot3( Y_helio_after(:,1)/AU, Y_helio_after(:,2)/AU, Y_helio_after(:,3)/AU, '-','LineWidth',2 );
hold on

% legend_entries{end+1}=[''];
legend_entries{end+1}=['Δ_',num2str(i)];
end

legend(legend_entries);

%% ex 2


% constants
mu_S = astroConstants(4);
mu_E= astroConstants(13);
AU = astroConstants(2);
R_E=astroConstants(23);

V_minus = [31.5; 5.2; 0.0]; % km/s
V_plus = [36.0; 0.0; 0.0]; % km/s
r_E = [0; -1; 0]*AU; % AU
V_P = [1; 0; 0]*sqrt(mu_S/norm(r_E));

v_inf_minus=V_minus-V_P;
v_inf_plus=V_plus-V_P;
delta=acos(dot(v_inf_minus,v_inf_plus)/(norm(v_inf_minus)*norm(v_inf_plus)));



        
    
    ecc_minus=@(rp) 1+(rp*(norm(v_inf_minus)^2)/mu_E);
    delta_minus=@(rp) 2*asin(1/ecc_minus(rp));

    ecc_plus=@(rp) 1+(rp*(norm(v_inf_plus)^2)/mu_E);
    delta_plus=@(rp) 2*asin(1/ecc_plus(rp));

    delta_fun=@(rp) delta_minus(rp)/2 + delta_plus(rp)/2;
    fun=@(rp) delta-delta_fun(rp);
rp=fzero(fun,R_E+2000);
% b=[1000:9000]';
% plot(b,fun(b));
h_atm=100; %km Karman line
if rp<R_E+h_atm 
    error('rp obtained is not physically feasible')
end 

a_minus=rp/(1-ecc_minus);
a_plus=rp/(1-ecc_plus);
vp_minus=sqrt(2*mu_E*(1/rp-1/(2*a_minus)));
vp_plus=sqrt(2*mu_E*(1/rp-1/(2*a_plus)));
