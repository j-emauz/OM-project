function [a_repeating, k, m]= RepeatingGroundTracks (rev,omega_E,T, mu_E,parameter)
% RepeatingGroundTracks Computes the semi-major axis, mean motion, and number of rotations
% of a satellite in order to achieve repeating ground tracks
%
% INPUT:
% rev: number of rotations of the planet or revolutions of the satellite
% omega_E: angular velocity of the planet [rad/s]
% T: orbital period of the satellite [s]
% mu_E: gravitational constant of the planet [km^3/s^2]
% parameter: 0 if the rotations of the planet are known, 1 if the revolutions of the satellite are known
%
% OUTPUT:
% a_repeating: semi-major axis of the orbit that achieves repeating ground tracks [km]
% k: number of rotations of the planet
% m: number of revolutions of the satellite
%
% Author:
% Name: Mariangela Testa, Oleksii Stepaniuk, João Emauz, Saverio Franzese
% Email: mariangela.testa@mail.polimi.it, oleksii.stepaniuk@mail.polimi.it,
% joao.emauz@mail.polimi.it, saverio.franzese@mail.polimi.it
n_SC=2*pi/T;

% Compute k from m or viceversa:
if parameter==0
     m=rev;
     k=m*n_SC/omega_E;
end 

if parameter==1 
    k=rev;
    m=k*omega_E/n_SC;
end 

% Repeating ground tracks can be obtained by choosing 
% an orbit with a mean motion such that:
n=omega_E*k/m;

%From which:

a_repeating=(mu_E/n^2)^(1/3);

end 
    