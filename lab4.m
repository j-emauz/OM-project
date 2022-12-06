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

up = [0; 0; 1]; %under planet
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

V_P = [0; 1; 0]*sqrt(mu_E/r_E);
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

% behind
r0=[-5*R_E; -Delta; 0 ];
y0=[r0; v_inf_minus];
% T=2*pi*sqrt( a^3/mu_E ); % Orbital period [1/s]
tspan = linspace( 0, T,1000);
     % Set options for the ODE solver
options = odeset( 'RelTol', 1e-14, 'AbsTol', 1e-14 );
[ Tb, Yb ] = ode113( @(t,y) ode_2bp(t,y,mu_E), tspan, y0, options );

% in front
r0=[-5*R_E; Delta; 0 ];
y0=[r0;v_inf_minus];
% T=2*pi*sqrt( a^3/mu_E ); % Orbital period [1/s]
tspan = linspace( 0, T,1000);
     % Set options for the ODE solver
options = odeset( 'RelTol', 1e-14, 'AbsTol', 1e-14 );
[ Tf, Yf ] = ode113( @(t,y) ode_2bp(t,y,mu_E), tspan, y0, options );

%% 3d plotting
figure(1)
earth_sphere;
grid on
hold on;
plot3( Yu(:,1), Yu(:,2), Yu(:,3), 'r-','LineWidth',2);

plot3( Yb(:,1), Yb(:,2), Yb(:,3), 'g-', 'LineWidth',2);

plot3( Yf(:,1), Yf(:,2), Yf(:,3), 'b-','LineWidth',2 );
legend('Flyby under the planet', 'Flyby behind the planet','Flyby in front the planet')