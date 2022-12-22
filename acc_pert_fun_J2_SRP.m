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
     [rr, vv] = par2car(s(1), s(2), s(3), s(4), s(5), s(6), mu_S);
%     
     r_dir = rr./norm(rr);
     w_dir = cross(rr,vv)./norm(cross(rr,vv));
     s_dir = cross(w_dir, r_dir);

     R_rsw2eci = [r_dir(:), w_dir(:), s_dir(:)];
    
    tilt = deg2rad(23.45);
    R_eci2sce = [1   0   0;                %ECI -> Sun-Centered-Ecliptic
          0  cos(tilt)   sin(tilt);
         0  -sin(tilt)  cos(tilt)];
   %R_sce2eci1 = [1, 0, 0];
   %R_sce2eci2 = [0 cos(tilt) -sin(tilt)];
   %R_sce2eci3 = [0 sin(tilt) cos(tilt)];

    r_sc_Sun= - R_eci2sce'*(r_S_E); 
     %r_sc_Sun1 = R_sce2eci1*(r_E_S);
     %r_sc_Sun2 = R_sce2eci2*(r_E_S);
     %r_sc_Sun3 = R_sce2eci3*(r_E_S);
     %r_sc_Sun = [r_sc_Sun1; r_sc_Sun2; r_sc_Sun3]
     %r_sc_Sun = r_E_S;

    p = s(1)*(1-s(2)^2);
    r = p/(1 + s(2)*cos(s(6)));
    v = sqrt(2*mu/r - mu/s(1));
    h = sqrt(p*mu);
    
    acc_pert_vec_J2 = -3/2*J2*mu*R^2/r^4.*[1-3*(sin(s(3)))^2*sin((s(6)+s(5)))^2; (sin(s(3)))^2*sin(2*(s(6)+s(5))); sin(2*s(3))*sin((s(6)+s(5)))];
    acc_pert_SRP= P_sun*(AU^2)/(norm(r_sc_Sun)^2)*Cr*AMR*10^-3;
    acc_pert_vec_SRP=-acc_pert_SRP*(r_sc_Sun./norm(r_sc_Sun));
    
    acc_pert_vec_SRP = R_rsw2eci'*acc_pert_vec_SRP;
    acc_pert_vec=acc_pert_vec_J2 + acc_pert_vec_SRP;

    % A_rot=h/(p*v).*[s(2)*sin(s(6)),  - (1+s(2)*cos(s(6)));
    %                    1+s(2)*cos(s(6)),   s(2)*sin(s(6))];
    % a_tn = A_rot\[a_rsw(1); a_rsw(2)];
    % acc_pert_vec=[a_tn; a_rsw(3)];
end