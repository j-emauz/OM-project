function [dv_tot,VF,car_2,t1,dt] = dv_arc1(t1, t2, p1, p2, mu)

dt = t2 - t1;
dt = dt*24*3600;
[kep_1,~] = uplanet(t1, p1);
[kep_2,~] = uplanet(t2, p2);

[car_1, v_1] = par2car(kep_1(1), kep_1(2), kep_1(3), kep_1(4), kep_1(5), kep_1(6), mu);
[car_2, v_2] = par2car(kep_2(1), kep_2(2), kep_2(3), kep_2(4), kep_2(5), kep_2(6), mu);

[A,P,E,ERROR,VI,VF,TPAR,THETA] = lambertMR(car_1, car_2, dt, mu, 0, 0, 0 );

dv_tot =norm(VI.' - v_1);

end