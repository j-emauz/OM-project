function [rp, Delta_vp,vp_minus,vp_plus,v_inf_minus,v_inf_plus] = PGA (V_P,V_minus,V_plus,rp_min,mu)
% PGA: Calculates the optimal planetary flyby radius and corresponding velocities.
%
% INPUT:
% V_P = Velocity of the planet [km/s]
% V_minus = Initial velocity of spacecraft [km/s]
% V_plus = Final velocity of spacecraft [km/s]
% rp_min = Minimum flyby radius [km]
% mu = Gravitational parameter of the central body [km^3/s^2]
%
% OUTPUT:
% rp = Optimal planetary flyby radius [km]
% Delta_vp = Powered flyby delta velocity [km/s]
% vp_minus = Velocity of spacecraft before flyby [km/s]
% vp_plus = Velocity of spacecraft after flyby [km/s]
% v_inf_minus = Flyby asymptotic velocity in the direction of V_minus [km/s]
% v_inf_plus = Flyby asymptotic velocity in the direction of V_plus [km/s]
%
% USAGE:
% [rp, Delta_vp,vp_minus,vp_plus,v_inf_minus,v_inf_plus] = PGA(V_P,V_minus,V_plus,rp_min,mu)
%
% Author:
% Name: Mariangela Testa, Oleksii Stepaniuk, João Emauz, Saverio Franzese
% Email: mariangela.testa@mail.polimi.it, oleksii.stepaniuk@mail.polimi.it,
% joao.emauz@mail.polimi.it, saverio.franzese@mail.polimi.it
options = optimset('TolFun',1e-14,'Display','off');


v_inf_minus=V_minus-V_P;
v_inf_plus=V_plus-V_P;
delta=acos(dot(v_inf_minus,v_inf_plus)/(norm(v_inf_minus)*norm(v_inf_plus)));
          
    ecc_minus=@(rp) 1+(rp*(norm(v_inf_minus)^2)/mu);
    delta_minus=@(rp) 2*asin(1/ecc_minus(rp));

    ecc_plus=@(rp) 1+(rp*(norm(v_inf_plus)^2)/mu);
    delta_plus=@(rp) 2*asin(1/ecc_plus(rp));

    delta_fun=@(rp) delta_minus(rp)/2 + delta_plus(rp)/2;
    fun=@(rp) delta-delta_fun(rp);
    
rp=fzero(fun,[rp_min, 10^10], options);

a_minus=rp/(1-ecc_minus(rp));
a_plus=rp/(1-ecc_plus(rp));
vp_minus=sqrt(2*mu*(1/rp-1/(2*a_minus)));
vp_plus=sqrt(2*mu*(1/rp-1/(2*a_plus)));
Delta_vp=abs(vp_plus-vp_minus);

if rp < rp_min || isnan(rp)
    Delta_vp=NaN;
end 

end