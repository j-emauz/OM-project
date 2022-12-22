function acc_pert_vec = acc_pert_fun_J2_SRP( t, s, mu,J2,Cr,AU,R,r_sc_Sun,v_SC)
 % Evaluate the perturbing accelerations in a given
    % RSW


    

    p = s(1)*(1-s(2)^2);
    r = p/(1 + s(2)*cos(s(6)));
    v = sqrt(2*mu/r - mu/s(1));
    h = sqrt(p*mu);
    
    acc_pert_vec_J2 = -3/2*J2*mu*R^2/r^4.*[1-3*(sin(s(3)))^2*sin((s(6)+s(5)))^2; (sin(s(3)))^2*sin(2*(s(6)+s(5))); sin(2*s(3))*sin((s(6)+s(5)))];
    acc_pert_SRP= P_sun*(AU^2)/(vecnorm(r_sc_Sun,2,2)^2)*Cr*AMR;
    acc_pert_vec_SRP=-acc_pert_vec_SRP*(r_sc_Sun./vecnorm(r_sc_Sun,2,2));
    
    acc_pert_vec_SRP
    acc_pert_vec=acc_pert_vec_J2+acc_pert_vec_SRP;

    % A_rot=h/(p*v).*[s(2)*sin(s(6)),  - (1+s(2)*cos(s(6)));
    %                    1+s(2)*cos(s(6)),   s(2)*sin(s(6))];
    % a_tn = A_rot\[a_rsw(1); a_rsw(2)];
    % acc_pert_vec=[a_tn; a_rsw(3)];
end