function acc_pert_vec = acc_pert_fun( t, s, mu,J2,R)
% Evaluate the perturbing accelerations in a given
% TNH

p = s(1)*(1-s(2)^2);
r = p/(1 + s(2)*cos(s(6)));
v = sqrt(2*mu/r - mu/s(1));
h = sqrt(p*mu);

a_rsw = -3/2*J2*mu*R^2/r^4*[1-3*(sin(s(3)))^2*sin((s(6)+s(5)))^2; (sin(s(3)))^2*sin(2*(s(6)+s(5))); sin(2*s(3))*sin((s(6)+s(5)))];

A_rot=h/(p*v).*[s(2)*sin(s(6)),  - (1+s(2)*cos(s(6)));
                   1+s(2)*cos(s(6)),   s(2)*sin(s(6))];
a_tn = A_rot\[a_rsw(1); a_rsw(2)];
acc_pert_vec=[a_tn; a_rsw(3)];
end