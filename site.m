%% SITE
%%
function R_site = site( sitlat,sitalt,lst )
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%-
%
%  Use         : R_site = site( sitlat,sitalt,lst )
%
%  The site function computes the Earth Centered Inertial (IJK) position
%  vector of a site on the Earth.
%
%  Author      : C2C Alekos Michael  USAFA  719-333-4529  16 Feb 2016
%
%  Inputs      :
%    sitlat      - Site Latitude                               (rad)
%    sitalt      - Site Altitude                                (km) 
%    lst         - Local Sidereal Time                    0.0 to 2Pi rad
%
%  Outputs     :
%    R_site      - Site Earth Centered Inertial position vector (km)
%
%  Locals      :
%    X           - Projection of R_site in x direction      
%    Z           - Projection of R_site in z direction      
%
%  Globals     :
%    Rad         - Conversion from degrees to radians
%    MU          - Earth Gravitational Constant           (km^3/s^2)
%    OmegaEarth  - Angular Speed of Earth                    (rad/s)
%    EEarth      - Eccentricity of Earth               
%    RE          - Radius of Earth                              (km)
%
%  Coupling    :
%    None.
%
%  References  :
%    None.
%
%  Documenation:
%    None.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%-
%Initialize Program Setup and Define Globals/Constants
wgs84data
global MU
global OmegaEarth
global Rad
global EEarth
global RE

%Calculate Projection of R_site in x direction
X=((RE/(sqrt(1-EEarth^2*sin(sitlat)^2)))+sitalt)*cos(sitlat);
%Calculate Projection of R_site in z direction
Z=(((RE*(1-EEarth^2))/(sqrt(1-EEarth^2*sin(sitlat)^2)))+sitalt)*sin(sitlat);
%Calculate R_site
R_site=[X*cos(lst),X*sin(lst),Z];

%End function
end

