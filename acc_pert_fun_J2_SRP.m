function acc_pert_vec = acc_pert_fun_J2_SRP( t, s, mu,J2,R,  initial_date, AMR, Cr, Perturbations)
% INPUT:
% t - time (seconds)
% s - state vector at current time (keplerian elements [a,e,i,Om,om,th]')
% mu - gravitational parameter of central body
% J2 - zonal coefficient of the central body
% R - radius of central body
% initial_date - vector of [yyyy, mm, dd, hh, mm, ss] representing initial time
% AMR - solar to mass ratio used to calculate solar radiation pressure
% Cr - reflectivity
% Perturbations - flag indicating which perturbations to include: 0: J2, 1: SRP, 2: J2+SRP
%
% OUTPUT:
% acc_pert_vec - acceleration vector due to perturbations in RSW (radial
% -transversal-out of plane) frame
%
% Functions used:
% astroConstants - returns astronomical constants
% date2mjd2000 - converts date to modified Julian date
% uplanet - returns position and velocity of planet in Solar System
% par2car - converts keplerian elements to Cartesian elements
    AU = astroConstants(2);% Astronomical Unit
    mu_S = astroConstants(4); % Gravitational parameter of the Sun
    P_sun = 4.57*10^-6; % Solar radiation pressure at 1 AU
    R_E = astroConstants(23); %radius of earth

    % Calculate time in days since initial_date
     initial_time=date2mjd2000(initial_date);
     t_days=t/3600/24; %days since initial time

     % Get position of Earth in Sun centered ecliptic frame
     [kep_E,~]=uplanet(initial_time+t_days,3);
     [r_S_E,~] = par2car(kep_E(1),kep_E(2),kep_E(3),kep_E(4),kep_E(5),kep_E(6),mu_S);

     % Get position of spacecraft in ECI frame
     [r_E_sc,~] = par2car(s(1), s(2), s(3), s(4), s(5), s(6), mu);
     
     % Calculate rotation matrices for transformation to RSW frame
    R_1_i=[1 0 0; 0 cos(s(3)) sin(s(3)); 0 -sin(s(3)) cos(s(3))];
    R_3_Om=[cos(s(4)) sin(s(4)) 0; -sin(s(4)) cos(s(4)) 0; 0 0 1];
    R_3_om_th=[cos(s(5)+s(6)) sin(s(5)+s(6)) 0; -sin(s(5)+s(6)) cos(s(5)+s(6)) 0; 0 0 1];
    
    % Calculate tilt of Earth's equator relative to the celestial equator
    tilt = deg2rad(23.45);
    % Calculate rotation matrix to transform from ECI to Sun-Centered-Ecliptic frame
    R_eci2sce = [1   0   0;              
          0  cos(tilt)   sin(tilt);
         0  -sin(tilt)  cos(tilt)];

    % Calculate vector of spacecraft to sun in ECI frame
    r_sc_Sun= - R_eci2sce'*(r_S_E); 
    
    % Calculate semi-latus rectum and r of spacecraft orbit
    p = s(1)*(1-s(2)^2);
    r = p/(1 + s(2)*cos(s(6)));
    

    % Calculate acceleration due to J2 perturbation in RSW frame
    acc_pert_vec_J2 = -3/2*J2*mu*R^2/r^4.*[1-3*(sin(s(3)))^2*sin((s(6)+s(5)))^2; (sin(s(3)))^2*sin(2*(s(6)+s(5))); sin(2*s(3))*sin((s(6)+s(5)))];
    % Calculate acceleration due to SRP perturbation in cartesian
    % coordinates
    acc_pert_SRP= P_sun*(AU^2)/(norm(r_sc_Sun)^2)*Cr*AMR*10^-3;
    acc_pert_vec_SRP=-acc_pert_SRP*(r_sc_Sun./norm(r_sc_Sun));

    % Transform acceleration due to SRP perturbation to RSW frame
    acc_pert_vec_SRP = R_3_om_th*R_1_i*R_3_Om*acc_pert_vec_SRP;

    % If spacecraft is in Earth's shadow, set acceleration due to SRP to zero
    r_E_sc_norm = norm(r_E_sc);
    r_E_sun_norm = norm(r_sc_Sun);
    theta_tot = acos(dot(r_E_sc,r_sc_Sun)/(r_E_sc_norm*r_E_sun_norm));
    theta_1 = acos(R_E/r_E_sc_norm);
    theta_2 = acos(R_E/r_E_sun_norm);
    if(theta_1+theta_2<=theta_tot)
         a_pert_vec_SRP = zeros(3,1);
     end
     %if sign(dot(r_E_sc./norm(r_E_sc),-r_sc_Sun./norm(r_sc_Sun)))==1 && (acos(dot(r_E_sc./norm(r_E_sc),-r_sc_Sun./norm(-r_sc_Sun)))<=70*pi/180 || acos(dot(r_E_sc./norm(r_E_sc),-r_sc_Sun./norm(-r_sc_Sun)))>= -70*pi/180)
     %   acc_pert_vec_SRP=zeros(3,1);
     %end
     
     % Select which perturbation to consider based on Perturbations flag
     if Perturbations==0
         acc_pert_vec=acc_pert_vec_J2;
     elseif Perturbations == 1
         acc_pert_vec=acc_pert_vec_SRP;
     else 
         acc_pert_vec=acc_pert_vec_J2 + acc_pert_vec_SRP;
     end

end