% Set your initial time and state (Cartesian or Keplerian): t0, s0
% Set your integration time span: tspan
% Set physical parameters (e.g., gravitational parameter of the primary, J2, etc.)
% Set ODE solver options (e.g., RelTol, AbsTol): options
% Numerical integration of the equations of motion
clc
close all
clear

mu_E=astroConstants(13);
R_E=astroConstants(23);
J2=astroConstants(9);
kep0 = [7571, 0.01, 87.9 *pi/180, 180 *pi/180, 180 *pi/180, 0 *pi/180];
[car0, v] = par2car(kep0(1), kep0(2), kep0(3), kep0(4), kep0(5), kep0(6), mu_E);
T = 2*pi*sqrt( kep0(1)^3/mu_E ); % Orbital period [1/s]
tspan = linspace( 0, 100*T, 10001 );
options = odeset( 'RelTol', 1e-13, 'AbsTol', 1e-14 );
[ T_Gauss, S_Gauss ] = ode45(@(t,s) eq_motion(t,s, @(t,s) acc_pert_fun(t,s,mu_E,J2,R_E), mu_E), tspan, kep0, options );
% Analyse and plot the results


options = odeset( 'RelTol', 1e-13, 'AbsTol', 1e-14 );
[ T_J2, S_J2 ] = ode113(@(t,s) perturbed_ode_2bp(t,s, mu_E, J2, R_E), tspan, [car0;v], options );
for k=1:size(S_J2,1)
   [a(k),e(k),i(k),OM(k),om(k),th(k),~]=car2par(S_J2(k,1:3)', S_J2(k,4:6)',mu_E);
%[a,e,i,OM,om,th, ee]=car2par(rr,vv,mu)
%@(t_tt, y_tt) ode_2bp(t_tt, y_tt, mu), tspan, y0_t, options);
end



%% plotting
hold on
N = 100;
amean = movmean(S_Gauss(:, 1), N);
plot(T_Gauss/T, amean);
plot(T_J2/T,a);
plot(T_Gauss/T,S_Gauss(:,1));
legend('average','Cartesian','Gauss')

figure
plot(T_Gauss/T,a-S_Gauss(:,1));
%% 
figure
N = 500;
emean = movmean(S_Gauss(:, 2), N);
plot(T_Gauss/T,S_Gauss(:,2));
hold on
plot(T_J2/T,e);
plot(T_Gauss/T, emean);
legend('Gauss','Cartesian','average')
%%
figure
N = 500;
emean = movmean(S_Gauss(:, 2), N);
plot(T_Gauss/T,S_Gauss(:,2));
hold on
plot(T_J2/T,e);
plot(T_Gauss/T, emean);
legend('Gauss','Cartesian','average')