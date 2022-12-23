function dy = perturbed_ode_2bp_SRP( t, y, mu, J2, R, initial_date, AMR, Cr, Perturbations)
%perturbations 0 - J2, 1- SRP, 2-SRP+J2

% Position and velocity
r = y(1:3);
v = y(4:6);
% Distance from the primary
rnorm = norm(r);
% Set the derivatives of the state

AU = astroConstants(2);
mu_S = astroConstants(4);
P_sun = 4.57*10^-6; %in N/m^2


 initial_time=date2mjd2000(initial_date);
 t_days=t/3600/24; %days since initial time
 [kep_E,~]=uplanet(initial_time+t_days,3);
 [r_S_E,~] = par2car(kep_E(1),kep_E(2),kep_E(3),kep_E(4),kep_E(5),kep_E(6),mu_S);

    
  tilt = deg2rad(23.45);
  R_eci2sce = [1   0   0;                %ECI -> Sun-Centered-Ecliptic
        0  cos(tilt)   sin(tilt);
        0  -sin(tilt)  cos(tilt)];
  r_sc_Sun= - R_eci2sce'*(r_S_E); 
  % sign
  acc_pert_SRP= P_sun*(AU^2)/(norm(r_sc_Sun)^2)*Cr*AMR*1e-3;
  acc_pert_vec_SRP= - acc_pert_SRP*(r_sc_Sun./norm(r_sc_Sun));
r_E_sc=r;
       if sign(dot(r_E_sc./norm(r_E_sc),-r_sc_Sun./norm(r_sc_Sun)))==1 && (acos(dot(r_E_sc./norm(r_E_sc),-r_sc_Sun./norm(-r_sc_Sun)))<=70*pi/180 || acos(dot(r_E_sc./norm(r_E_sc),-r_sc_Sun./norm(-r_sc_Sun)))>= -70*pi/180)
             % if the two vectors are in the same direction, within a cone of 160 deg,
             % the SC is in the shadow cone
            acc_pert_vec_SRP=zeros(3,1);
       end
       
       


a=3/2*J2*mu*R^2/rnorm^4*[r(1)/rnorm*(5*r(3)^2/rnorm^2-1); r(2)/rnorm*(5*r(3)^2/rnorm^2-1); r(3)/rnorm*(5*r(3)^2/rnorm^2-3)   ];

if(Perturbations == 0)
    acc_pert_vec_SRP = 0;
elseif Perturbations == 1
    a = 0;
end

dy = [ v
(-mu/rnorm^3)*r + a + acc_pert_vec_SRP];


end
