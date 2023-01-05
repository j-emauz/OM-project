function [a,e,i,OM,om,th, ee]=car2par(rr,vv,mu)
% car2par: Convert ECI Position and Velocity Vectors to Keplerian Orbital Elements
%
% This function converts the position and velocity vectors of an object in an
% Earth-centered inertial reference frame (ECI) to Keplerian orbital elements.
%
% The input parameters are the position vector "rr" and velocity vector "vv" in
% the ECI frame, and the gravitational constant "mu" for the attractor body.
% The output parameters are the semi-major axis "a", eccentricity "e",
% inclination "i", right ascension of the ascending node "OM", argument of
% perigee "om", and true anomaly "th" in radians. The function also prints the
% angles in degrees to the screen. If the orbit is equatorial, the right
% ascension of the ascending node is set to zero and the argument of perigee is
% set to the angle between the x-axis and the eccentricity vector in the xy-plane.
% If the orbit is circular, the argument of perigee and true anomaly are set to
% "NaN".
%
% Author:
% Name: Mariangela Testa
% Email: mariangela.testa@mail.polimi.it
%
% Last review date: January 5, 2023
%
% Usage: [a,e,i,OM,om,th, ee]=car2par(rr,vv,mu)
%
% Inputs:
% rr: COLUMN vector r in X,Y,Z components [km]
% vv: COLUMN vector v in X,Y,Z components [km/s]
% mu: planetary gravitational constant G*m1, with m1 mass of the attractor
% for the Earth: mu=398600 [km^3/s^2]
%
% Outputs:
% KEPLERIAN PARAMETERS
% a: semi-major axis
% e: eccentricity
% i: inclination of orbital plane
% OM: right ascension of ascending node (RAAN)
% om: argument of perigee
% th: true anomaly
%

%!!! NB: vectors are indicated with double letters, magnitudes with single letters
% eg. rr and r
%!!! NB: om and th are NaN in the case of circular orbits

%INPUT
% rr: COLUMN vector r in X,Y,Z components [km]
% vv: COLUMN vector v in X,Y,Z components [km/s]
% mu: planetary gravitational constant G*m1, with m1 mass of the attractor
% for the Earth: mu=398600 [km^3/s^2]

%OUTPUT
% KEPLERIAN PARAMETERS
% a: semi-major axis
% e: eccentricity
% i: inclination of orbital plane
% OM: right ascension of ascending node (RAAN)
% om: argument of perigee
% th: true anomaly

xx=[1,0,0]';
yy=[0,1,0]';
zz=[0,0,1]';

% 1): norms --------------------------------------------------------------
r=norm(rr);
v=norm(vv);

% 2): a ------------------------------------------------------------------
a=1/(2/r-v^2/mu); %[km]

% 3): hh vector, h magnitude -----------------------------------------------
hh=cross(rr,vv);
h=norm(hh);

% 4): ee vector, e magnitude -----------------------------------------------
ee=1/mu*cross(vv,hh)-rr./r;
e=norm(ee);

% 5): n ascending node unit vector --------------------------------------
nn=cross(zz,hh)/norm(cross(zz,hh));

% 6): i inclination angle ----------------------------------------------
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



% 10): CHECKS for particular orbits

if i==0 || i==pi %% identify equatorial orbits
    disp('Orbita equatoriale');
    OM=0;
    om=acos((xx'*ee)./e);
    if ee'*yy < 0
    om=2*pi-om;
    end
end

if e<=1e-10 %% identify circular orbits
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



    