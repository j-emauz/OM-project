function [rr, vv] = par2car(a, e, i, OM, om, theta, mu)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% DESCRIPTION par2car.m:                                                                     
% Converte le coordinate kepleriane classiche dell'orbita d'interesse nei           
% vettori r e v scritti rispetto sistema di riferimento celeste (ECI).                                                                                       
%                                                                                                                                   
% INPUT:                                                                           
% - a [km]:          semiasse maggiore dell'ellisse orbitale                           
% - e:               eccentricità dell'orbita                                          
% - i [rad]:         inclinazione dell'orbita. Corrisponde all'angolo 
%                    compreso tra l'asse Z del sistema ECI e la normale 
%                    al piano orbitale        
% - OM [rad]:        ascensione retta del nodo ascendente. Corrisponde 
%                    all'angolo tra il versore N che indica la direzione 
%                    del nodo ascendente   
%                    (intersezione orbita - piano equatoriale) e l'asse X 
%                    del sistema ECI (direzione dell'equinozio di 
%                    primavera)            
% - om [rad]:        anomalia del pericentro. Corrisponde all'angolo 
%                    compreso tra la direzione del pericentro e il versore 
%                    N descritto prima     
% - theta [rad]:     anomalia vera. Angolo che mi da la posizione 
%                    del satellite sull'orbita
% - mu [km^3 / s^2]: prodotto di G (costante di gravitazione universale) e
%                    la massa dell'attrattore principale
%
% OUTPUT:                                                                           
% - rr [km]:         vettore che mi da la posizione del satellite nello 
%                    spazio rispetto al sistema ECI                                       
% - vv [km/s]:       vettore che mi da la velocità del satellite nel punto 
%                    indicato da rr rispetto al sistema ECI                                  
% Author:
% Name: Mariangela Testa
% Email: mariangela.testa@mail.polimi.it
%                                                                     
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% 1 Calcolo semilato retto

p = a * (1 - e^2); % [km]                                                  % semilato retto

%% 2 Calcolo modulo del vettore r

r = p / (1 + e*cos(theta)); % [km]                                         % modulo del vettore r

%% 3 Calcolo vettore r in sistema perifocale

rr_pf = r * [cos(theta); sin(theta); 0]; % [km]                            % vettore r scritto nel sistema pf

%% 4 Calcolo del vettore v in sistema perifocale

vr = sqrt(mu/p) * e * sin(theta); % [km / s]                               % modulo del vettore velocità radiale
vt = sqrt(mu/p) * (1 + e*cos(theta)); % [km / s]                           % modulo del vettore velocità tangenziale

iir = [cos(theta); sin(theta); 0];                                         % versore radiale nel sistema pf
iit = [-sin(theta); cos(theta); 0];                                        % versore tangente nel sistema pf

vv_pf = vr * iir + vt * iit; % [km / s]                                    % vettore velocità nel sistema pf

%% 5 Definizione matrici di rotazione

R3_OM = [cos(OM), sin(OM), 0; -sin(OM), cos(OM), 0; 0, 0, 1];               
R1_i = [1, 0, 0; 0, cos(i), sin(i); 0, -sin(i), cos(i)];                   
R3_om = [cos(om), sin(om), 0; -sin(om), cos(om), 0; 0, 0, 1];              

%% 6 Matrice di rotazione complessiva

R = R3_om * R1_i * R3_OM;                                                   

%% 7 Calcolo vettori v e r nel sistema ECI

rr = R' * rr_pf; % [km]
vv = R' * vv_pf; % [km / s]