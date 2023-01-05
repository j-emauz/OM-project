function ds = eq_motion( t, s, acc_pert_fun_J2_SRP, mu)
% ODE system for perturbed two-body problem
    % INPUTS:
    % - t: scalar, time
    % - s: 6x1 vector, state vector (keplerian elements [a,e,i,Om,om,th]') at time t 
    % - acc_pert_fun_J2_SRP: function handle, function that calculates perturbing accelerations
    % - mu: scalar, gravitational parameter of central body
    % OUTPUTS:
    % - ds: 6x1 vector, time derivative of state vector at time t
    % FUNCTIONS USED:
    % - acc_pert_fun_J2_SRP: function handle, calculates perturbing accelerations
    
% Authors
% Name: Mariangela Testa, Oleksii Stepaniuk, João Emauz, Saverio Franzese
% Email: mariangela.testa@mail.polimi.it, oleksii.stepaniuk@mail.polimi.it,
% joao.emauz@mail.polimi.it, saverio.franzese@mail.polimi.it

    % Evaluate the perturbing accelerations
    acc_pert = acc_pert_fun_J2_SRP( t, s );

    
    p = s(1)*(1-s(2)^2);
    r = p/(1 + s(2)*cos(s(6)));
    h = sqrt(p*mu);
    
    % Calculate time derivatives of keplerian elements
    da=2*(s(1))^2/h*(s(2)*sin(s(6))*acc_pert(1)+p/r*acc_pert(2));
    de=1/h*(p*sin(s(6))*acc_pert(1)+acc_pert(2)*((p+r)*cos(s(6))+r*s(2)));
    di=r/h*cos(s(6)+s(5))*acc_pert(3);
    dOm=r/(h*sin(s(3)))*sin(s(6)+s(5))*acc_pert(3);
    dom=1/(h*s(2))*(-p*cos(s(6))*acc_pert(1)+(p+r)*sin(s(6))*acc_pert(2))-r*sin(s(6)+s(5))*cos(s(3))*acc_pert(3)/(h*sin(s(3)));
    df=h/r^2+1/(s(2)*h)*(p*cos(s(6))*acc_pert(1)-(p+r)*sin(s(6))*acc_pert(2));
    
    ds = [da;de;di;dOm;dom;df];


end