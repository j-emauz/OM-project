% 1.a
clear
close all

% constants
mu_S = astroConstants(4);
mu_E= astroConstants(13);
AU = astroConstants(2);

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