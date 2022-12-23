function acc_pert_vec = acc_pert_fun_J2_SRP( t, s, mu,J2,R)
 % Evaluate the perturbing accelerations in a given
    % RSW
    Cr = 1;
    AU = astroConstants(2);
    AMR = 10;
    mu_S = astroConstants(4);
    P_sun = 4.5*10^-6;

     initial_date=[2022,03,21,12,0,0];
     initial_time=date2mjd2000(initial_date);
     t_days=t/3600/24; %days since initial time
     [kep_E,~]=uplanet(initial_time+t_days,3);
     [r_S_E,~] = par2car(kep_E(1),kep_E(2),kep_E(3),kep_E(4),kep_E(5),kep_E(6),mu_S);
     [r_E_sc,~] = par2car(s(1), s(2), s(3), s(4), s(5), s(6), mu);
     
     %r_dir = rr./norm(rr);
     %w_dir = cross(rr,vv)./norm(cross(rr,vv));
     %s_dir = cross(w_dir, r_dir);

    R_1_i=[1 0 0; 0 cos(s(3)) sin(s(3)); 0 -sin(s(3)) cos(s(3))];
    R_3_Om=[cos(s(4)) sin(s(4)) 0; -sin(s(4)) cos(s(4)) 0; 0 0 1];
    R_3_om_th=[cos(s(5)+s(6)) sin(s(5)+s(6)) 0; -sin(s(5)+s(6)) cos(s(5)+s(6)) 0; 0 0 1];

    

     %R_rsw2eci = [r_dir(:), s_dir(:), w_dir(:)];
    
    tilt = deg2rad(23.45);
    R_eci2sce = [1   0   0;                %ECI -> Sun-Centered-Ecliptic
          0  cos(tilt)   sin(tilt);
         0  -sin(tilt)  cos(tilt)];

    r_sc_Sun= - R_eci2sce'*(r_S_E); 
    % Earth - Sun vector in the ECI frame

    p = s(1)*(1-s(2)^2);
    r = p/(1 + s(2)*cos(s(6)));
    v = sqrt(2*mu/r - mu/s(1));
    h = sqrt(p*mu);
    
    acc_pert_vec_J2 = -3/2*J2*mu*R^2/r^4.*[1-3*(sin(s(3)))^2*sin((s(6)+s(5)))^2; (sin(s(3)))^2*sin(2*(s(6)+s(5))); sin(2*s(3))*sin((s(6)+s(5)))];
    acc_pert_SRP= P_sun*(AU^2)/(norm(r_sc_Sun)^2)*Cr*AMR*10^-3;
    acc_pert_vec_SRP=-acc_pert_SRP*(r_sc_Sun./norm(r_sc_Sun));
    
    %acc_pert_vec_SRP = R_rsw2eci'*acc_pert_vec_SRP;
    acc_pert_vec_SRP = R_3_om_th*R_1_i*R_3_Om*acc_pert_vec_SRP;

     if dot(r_E_sc./norm(r_E_sc),-r_sc_Sun./norm(r_sc_Sun))==1 && (acos(dot(r_E_sc./norm(r_E_sc),-r_sc_Sun./norm(-r_sc_Sun)))<=70*pi/180 || acos(dot(r_E_sc./norm(r_E_sc),-r_sc_Sun./norm(-r_sc_Sun)))>= -70*pi/180)
%              if the two vectors are in the same direction, within a cone of 160 deg,
%              the SC is in the shadow cone
        acc_pert_vec_SRP=zeros(3,1);
     end

     acc_pert_vec=acc_pert_vec_J2 + acc_pert_vec_SRP;

end