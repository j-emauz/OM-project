function dy = perturbed_ode_2bp_SRP( t, y, mu, J2, R, initial_date, AMR, Cr, Perturbations)
%INPUTS
% t: time
% y: state, y = [position; velocity]
% mu: gravitational parameter of the central body
% J2: zonal coefficient
% R: radius of the central body
% initial_date: initial date in the format [year, month, day, hour, minute, second]
% AMR: area-to-mass ratio of the spacecraft
% Cr: coefficient of reflectivity of the spacecraft
% Perturbations - flag indicating which perturbations to include: 0: J2, 1: SRP, 2: J2+SRP
%OUTPUTS
% dy: state derivative, dy = [velocity; acceleration]
% USED FUNCTIONS:
% - astroConstants: returns astronomical constants
% - date2mjd2000: converts a date to the modified Julian date (2000 epoch)
% - uplanet: computes the Keplerian elements of a planet in SCE frame
% - par2car: converts Keplerian elements to Cartesian coordinates

% Position and velocity
r = y(1:3);
v = y(4:6);
% Norm of the position vector
rnorm = norm(r);
% Set the derivatives of the state

AU = astroConstants(2); % Astronomical unit in meters
mu_S = astroConstants(4);  % Gravitational parameter of the sun
P_sun = 4.57*10^-6;% Solar radiation pressure at 1 AU

% Convert time to days since the initial date
 initial_time=date2mjd2000(initial_date);% Convert initial date to modified Julian date
 t_days=t/3600/24; %days since initial time

 % Get position of Earth in Sun centered ecliptic frame
 [kep_E,~]=uplanet(initial_time+t_days,3);
 [r_S_E,~] = par2car(kep_E(1),kep_E(2),kep_E(3),kep_E(4),kep_E(5),kep_E(6),mu_S);

  % Calculate tilt of Earth's equator relative to the celestial equator
  tilt = deg2rad(23.45);
  % Calculate rotation matrix to transform from ECI to Sun-Centered-Ecliptic frame
  R_eci2sce = [1   0   0;                %ECI -> Sun-Centered-Ecliptic
        0  cos(tilt)   sin(tilt);
        0  -sin(tilt)  cos(tilt)];
  % Calculate vector of spacecraft to sun in ECI frame   
  r_sc_Sun= - R_eci2sce'*(r_S_E); 

  % Compute the acceleration due to SRP
  acc_pert_SRP= P_sun*(AU^2)/(norm(r_sc_Sun)^2)*Cr*AMR*1e-3;
  acc_pert_vec_SRP= - acc_pert_SRP*(r_sc_Sun./norm(r_sc_Sun));
  r_E_sc=r;
  % If spacecraft is in Earth's shadow, set acceleration due to SRP to zero
  if sign(dot(r_E_sc./norm(r_E_sc),-r_sc_Sun./norm(r_sc_Sun)))==1 && (acos(dot(r_E_sc./norm(r_E_sc),-r_sc_Sun./norm(-r_sc_Sun)))<=70*pi/180 || acos(dot(r_E_sc./norm(r_E_sc),-r_sc_Sun./norm(-r_sc_Sun)))>= -70*pi/180)
       % if the two vectors are in the same direction, within a cone of 160 deg,
          % the SC is in the shadow cone
       acc_pert_vec_SRP=zeros(3,1);
  end
       
       

% Compute the acceleration due to J2
a=3/2*J2*mu*R^2/rnorm^4*[r(1)/rnorm*(5*r(3)^2/rnorm^2-1); r(2)/rnorm*(5*r(3)^2/rnorm^2-1); r(3)/rnorm*(5*r(3)^2/rnorm^2-3)   ];

% Select which perturbation to consider based on Perturbations flag
if(Perturbations == 0)
    acc_pert_vec_SRP = 0;
elseif Perturbations == 1
    a = 0;
end

% Compute the state derivative
dy = [ v
(-mu/rnorm^3)*r + a + acc_pert_vec_SRP];


end
