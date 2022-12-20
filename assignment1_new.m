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
rp_min=R_E*21 * 1.1;

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
    for j=1:length(tspan_GA)       
        [dv_1(i,j),V_SC_Saturn_1,V_Saturn, r_Saturn, t1(i,j),ToF1] = dv_arc1(tspan_dept(i), tspan_GA(j), p1, p2, mu_S);        
        % dv_1 manoeuver at Earth
        for k=1:length(tspan_arrt) 
            [kepNEO_arr,~,~] = ephNEO(tspan_arrt(k),86);
            [dv_2(i,j,k),V_SC_Saturn_2, t2(i,j,k),ToF2] = dv_arc2(tspan_GA(i), tspan_arrt(k), r_Saturn, kepNEO_arr, mu_S);
            [rp, Delta_vp(i,j,k)] = PGA (V_Saturn, V_SC_Saturn_1, V_SC_Saturn_2, rp_min,mu_Saturn);
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