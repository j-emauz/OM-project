function [a,e,i,OM,om,th, ee]=car2par(rr,vv,mu)
% [a,e,i,OM,om,th,ee]=car2par(rr,vv,mu)
% SCOPO: passo da sistema di rif. geocentrico inerziale (Earth-centered-inertial, ECI)
% dati i vettori r e v orientati nello spazio
% restituisce angoli in radianti, stampa a schermo gli angoli in gradi

%!!! NB: i vettori sono indicati con lettera doppia, i moduli singola
%es rr e r
%!!! NB: om e th sono NaN nel caso di orbite circolari


%INPUT
% rr: vettore r COLONNA in componenti X,Y,Z     [km]
% vv: vettore v COLONNA in componenti X,Y,Z     [km/s]
% mu: costante gravitazionale planetaria G*m1, con m1 massa dell'attrattore
%
% per la Terra: mu=398600 [km^3/s^2]

%OUTPUT
% PARAMETRI KEPLERIANI
% a: SEMIasse maggiore
% e: eccentricità
% i: inclinazione del piano orbitale
% OM: ascensione retta del nodo ascendente (RAAN)
% om: anomalia  del pericentro
% th: anomalia vera

xx=[1,0,0]';
yy=[0,1,0]';
zz=[0,0,1]';

% 1): norme --------------------------------------------------------------
r=norm(rr);
v=norm(vv);

% 2): a ------------------------------------------------------------------
a=1/(2/r-v^2/mu);  %[km]

% 3): hh vettore, h modulo -----------------------------------------------
hh=cross(rr,vv);
h=norm(hh);

% 4): ee vettore, e modulo -----------------------------------------------
ee=1/mu*cross(vv,hh)-rr./r;
e=norm(ee);

% 5): n versore del nodo ascendente --------------------------------------
nn=cross(zz,hh)/norm(cross(zz,hh));

% 6): i angolo inclinazione ----------------------------------------------
i=acos(hh'*zz/h); %acos(hz/h)

% 7): OM -----------------------------------------------------------------
OM=acos(xx'*nn);
if nn'*yy < 0
OM=2*pi-OM;
end

% 8): om -----------------------------------------------------------------
om=acos((nn'*ee)./e);
if ee'*zz < 0
om=2*pi-om;
end

% 9): th -----------------------------------------------------------------
th=acos((rr'*ee)/(r*e));
if rr'*vv < 0
    th=2*pi-th;
end



% 10): CONTROLLI per orbite particolari

if i==0 || i==pi %% individuo le orbite equatoriali
    disp('Orbita equatoriale');
    OM=0;
    om=acos((xx'*ee)./e);
    if ee'*yy < 0
    om=2*pi-om;
    end
end

if e<=1e-10 %% individuo le orbite circolari
    disp('Orbita circolare');
    om=NaN;
    th=NaN;
end


% % OUTPUT in deg
% disp(['i=',num2str(i*180/pi),'°']);
% disp(['OM=',num2str(OM*180/pi),'°']);
% disp(['om=',num2str(om*180/pi),'°']);
% disp(['th=',num2str(th*180/pi),'°']);


end



    