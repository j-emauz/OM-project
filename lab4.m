% 1.a
clear
close all

% constants
mu_S = astroConstants(4);
mu_E= astroConstants(13);
AU = astroConstants(2);

% data
v_inf_minus = [15.1; 0; 0]; % km/s
Delta = 9200; %km 
r_E= [1; 0; 0]*AU; % AU

up = [0; 0; 1]; %under planet
fp = [0; 1; 0]; %front planet
bp = [0; -1; 0]; %behind planet


uu = cross(Delta*up, v_inf_minus)/norm(cross(Delta*up, v_inf_minus));
uf = cross(Delta*fp, v_inf_minus)/norm(cross(Delta*fp, v_inf_minus));
ub = cross(Delta*bp, v_inf_minus)/norm(cross(Delta*bp, v_inf_minus));