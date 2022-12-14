function ds = eq_motion( t, s, acc_pert_fun, mu)
    % Evaluate the perturbing accelerations
    
    acc_pert = acc_pert_fun( t, s );
    % Evaluate the equations of motion (Cartesian or Keplerian),
    % function of t, s, acc_pert_vec, and parameters
    
    p = s(1)*(1-s(2)^2);
    r = p/(1 + s(2)*cos(s(6)));
    v = sqrt(2*mu/r - mu/s(1));
    h = sqrt(p*mu);
    
    
    % da = 2*s(1)^2*v/mu*acc_pert(1);
    % de = 1/v*(2*(s(2)+cos(s(6))*acc_pert(1)-r/s(1)*sin(s(6))*acc_pert(2)));
    % di = r*cos(s(6)+s(5))/h *acc_pert(3);
    % dOm = r*sin(s(6)+s(5))/(h*sin(s(3)))*acc_pert(3);
    % dom = 1/(s(2)*v)*(2*sin(s(6))*acc_pert(1) + (2*s(2)+r/s(1)*cos(s(6)))*acc_pert(2)) - r*sin(s(6)+s(5))*cos(s(3))/(h*sin(s(3)))*acc_pert(3);
    % df = h/r^2 - 1/(s(2)*v) * (2*sin(s(6))*acc_pert(1)+(2*s(2)+r/s(1)*cos(s(6)))*acc_pert(2));
    
    da=2*(s(1))^2/h*(s(2)*sin(s(6))*acc_pert(1)+p/r*acc_pert(2));
    de=1/h*(p*sin(s(6))*acc_pert(1)+acc_pert(2)*((p+r)*cos(s(6))+r*s(2)));
    di=r/h*cos(s(6)+s(5))*acc_pert(3);
    dOm=r/(h*sin(s(3)))*sin(s(6)+s(5))*acc_pert(3);
    dom=1/(h*s(2))*(-p*cos(s(6))*acc_pert(1)+(p+r)*sin(s(6))*acc_pert(2))-r*sin(s(6)+s(5))*cos(s(3))*acc_pert(3)/(h*sin(s(3)));
    df=h/r^2+1/(s(2)*h)*(p*cos(s(6))*acc_pert(1)-(p+r)*sin(s(6))*acc_pert(2));
    
    ds = [da;de;di;dOm;dom;df];


end