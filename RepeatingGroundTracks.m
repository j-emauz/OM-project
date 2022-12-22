function [a_repeating, k, m]= RepeatingGroundTracks (rev,omega_E,T, mu_E,parameter)
%GIVE AS IMPUT:
%Parameter=0 if the rotations of the planet are known;
%Parameter=1 if the revolutions of the satellite are known;

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
    