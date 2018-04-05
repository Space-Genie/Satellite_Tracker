%% Rvector
%%
function [ R_ijk ] = Rvector( n,ecc,inc,raan,argp,nu )
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%-
%
%  Rvector uses COEs to calculate the spacecraft IJK position vector.
%
%  Author      : C2C Alekos Michael  USAFA  719-333-4529  15 Apr 2016
%
%  Inputs      :
%   n               Mean Motion                                     rad/sec
%   ecc             Eccentricity                                   unitless
%   inc             Inclination                                     rad
%   raan            Right Ascension of Ascending Node               rad
%   argp            Argument of Periapsis                           rad
%   nu              True Anomaly                                    rad
%
%  Outputs     :
%   R_ijk           Spacecraft IJK Position Vector                  km
%
%  Locals      :
%   a               Semi-Major Axis                                 km
%   p               Semi-Latus Rectum                               km
%   R_pqw           Position Vector in PQW Frame                    km
%   
%  Globals     :
%   MU              Earth Gravitational Constant                   km^3/s^2
%
%  Coupling    :
%   axisrot         Performs rotation of specified angle about a desired
%                   axis
%
%  References  :
%   None.
%
%  Documenation:
%   None.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%-
% Define Globals
global MU

PQW=[cos(nu), sin(nu),0]; %Defines PQW vector
a=(MU/n^2)^(1/3); %updates semi-major axis
p=a*(1-ecc^2); %updates semi-latus rectum
R_pqw=(p/(1+ecc*cos(nu)))*PQW; %calculates position vector in PQW frame
temp1=axisrot(R_pqw,3,-argp); %first rotation
temp2=axisrot(temp1,1,-inc); %second rotation
R_ijk=axisrot(temp2,3,-raan); %final rotation and final answer for R vector
%in IJK frame

end

