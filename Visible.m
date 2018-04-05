%% Visible
%%
function [rho,az,el,vis,Alpha,Beta,X] = Visible( R_ijk,R_site,sitlat,lst,jd )
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%-
%
%  Visible determines if the observation location is in the dark and the
%  visibility criteria to include elevation > 10 deg, range < 1500 km, and
%  if the satellite is illuminated. 
%
%  Author      : C2C Alekos Michael  USAFA  719-333-4529  15 Apr 2016
%
%  Inputs      :
%   R_ijk           Spacecraft IJK Position Vector                  km
%   R_site          Site IJK Position Vector                        km
%   sitlat          Site Geodetic Latitude                          rad
%   lst             Site Local Sidereal Time                        rad
%   jd              Julian Date at Viewing Time			

%
%  Outputs     :
%   rho             Range                                           km
%   az              Azimuth                                         rad
%   el              Elevation                                       rad			
%   vis             Flag for Visibility                             0 or 1
%
%  Locals      :
%   Beta            Angle Between Rsun and R_site Vectors           rad
%   Alpha           Angle Between Rsun and R_ijk Vectors            rad
%   X               Distance of Satellite From Earth-Sun Axis       km
%   Rsun            Position Vector of the Sun in IJK Frame         km   
%
%  Globals     :
%   RE              Radius of Earth                                 km
%   Rad             Conversion from degrees to radians                    
%
%  Coupling    :
%   sun             Calculates the Geocentric Equatorial position vector
%                   for the Sun given the Julian Date. 
%   vecangle        Finds the angle between two vectors
%   mag             Computes magnitude of a vector
%   RhoAzEl         Finds range, azimuth, and elevation information
%
%  References  :
%   None.
%
%  Documenation:
%   None.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%-
% Define Globals
global RE Rad

[RSun, RtAsc, Decl] = sun (jd); %gets sun vector Rsun
Alpha = vecangle(RSun,R_site); %Finds angle between sun and site vectors
Beta = vecangle(RSun,R_ijk); %Finds angle between sun and ijk vectors
X=R_ijk*sin(Beta); %calculates distance from Sun-Earth axis
X=mag(X); %calculates magnitude of X distance
[rho,az,el] = RhoAzEl( R_ijk,R_site,sitlat,lst ); %Obtains rho, az, el
%information
vis=0; %indicates visible conditions as No or "0" unless proven otherwise
if (Alpha >= pi/2)&&(el >= 10*Rad)&&(rho <= 1500)&&(Beta <= pi/2)||(Alpha >= pi/2)&&(el >= 10*Rad)&&(rho <= 1500)&&(X >= RE)
        vis=1; %this if statement checks the visible criteria and sets 
        %visible conditions to Yes or "1" if the if statementsare true
end
end

