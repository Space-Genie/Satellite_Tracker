%% RhoAzEl
%%
function [ rho,az,el ] = RhoAzEl( R_ijk,R_site,sitlat,lst )
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%-
%
%  This procedure computes the topocentric range, azimuth, and elevation
%  from the site IJK vector and the satellite IJK position vector.
%
%  Author      : C2C Alekos Michael  USAFA  719-333-4529  15 Apr 2016
%
%  Inputs      :
%   R_ijk           Spacecraft IJK Position Vector                  km
%   R_site          Site IJK Position Vector                        km
%   sitlat          Site Geodetic Latitude                          rad
%   lst             Site Local Sidereal Time                        rad
%
%  Outputs     :
%   rho             Range                                           km
%   az              Azimuth                                         rad
%   el              Elevation                                       rad
%
%  Locals      :
%   colat           Colatitude                                      rad
%   rho_sez         Range Vector in SEZ frame                       km
%   rho_ijk         Range Vector in IJK frame                       km
%   rhoS            S component of rho                              km
%   rhoE            E component of rho                              km
%   rhoZ            Z component of rho                              km
%
%  Globals     :
%   RE              Radius of Earth                                 km
%   MU              Earth Gravitational Constant                   km^3/s^2
%
%  Coupling    :
%   axisrot         Performs rotation of specified angle about a desired
%                   axis
%   revcheck        Accomplishes a modulus function of x by modby.
%
%  References  :
%   None.
%
%  Documenation:
%   None.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%-
% Define Globals
global RE MU

colat=pi/2-sitlat; %calculates colat
R_site=R_site'; %rotates Rsite vector to match other vectors in code
rho_ijk=R_ijk-R_site; %finds rho IJK
temp1=axisrot(rho_ijk,3,lst); %first rotation
rho_sez=axisrot(temp1,2,colat); %second rotation and final value for RhoSEZ
rho=mag(rho_sez); %magnitude of rho vector
rhoS=rho_sez(1); %rho S component 
rhoE=rho_sez(2); %rho E component 
rhoZ=rho_sez(3); %rho Z component 
az=atan2(rhoE,-rhoS); %calculates az in correct quadrant
az=revcheck(az,2*pi); %ensures az is less than 2*pi
el=asin(rhoZ/rho); %calculates elevation

end

