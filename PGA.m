function [rp, Delta_vp] = PGA (V_P,V_minus,V_plus,rp_min,mu)

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
    
rp=fzero(fun,rp_min, options);

a_minus=rp/(1-ecc_minus(rp));
a_plus=rp/(1-ecc_plus(rp));
vp_minus=sqrt(2*mu*(1/rp-1/(2*a_minus)));
vp_plus=sqrt(2*mu*(1/rp-1/(2*a_plus)));
Delta_vp=vp_plus-vp_minus;

if rp < rp_min || isnan(rp)
    Delta_vp=NaN;
end 

end