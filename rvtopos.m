%% RVTOPOS
%%
function [ Rho_sez,Drho_sez ] = rvtopos( rho,az,el,drho,daz,del )
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%-
% 
%  Use         : [ Rho_sez,Drho_sez ] = rvtopos( rho,az,el,drho,daz,del )
% 
%  The rvtopos function converts a radar measurement to a range and 
%  range rate vector in the SEZ frame.
% 
%  Author      : C2C Alekos Michael  USAFA  719-333-4529  16 Feb 2016
% 
%  Inputs      :
%    rho         - Range                                        (km)
%    az          - Azimuth                                     (rad)
%    el          - Elevation                                   (rad)
%    drho        - Range Rate                                 (km/s)
%    daz         - Azimuth Rate                              (rad/s)
%    del         - Elevation Rate                            (rad/s)
% 
%  Outputs     :
%    Rho_sez     - Relative Position Vector in SEZ frame        (km)
%    Drho_sez    - Relative Velocity Vector in SEZ frame      (km/s)
% 
%  Locals      :
%    RhoS        - S component of Rho vector                    (km)
%    RhoE        - E component of Rho vector                    (km)
%    RhoZ        - Z component of Rho vector                    (km)
%    DrhoS       - S component of Drho vector                   (km)
%    DrhoE       - E component of Drho vector                   (km)
%    DrhoZ       - Z component of Drho vector                   (km)
%      
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

%Calculate S,E, and Z components of Rho_sez vector
RhoS=-rho*cos(el)*cos(az);
RhoE=rho*cos(el)*sin(az);
RhoZ=rho*sin(el);

%Calculate S,E, and Z components of Drho_sez vector
DrhoS=-drho*cos(el)*cos(az)+rho*del*sin(el)*cos(az)+rho*daz*cos(el)*sin(az);
DrhoE=drho*cos(el)*sin(az)-rho*del*sin(el)*sin(az)+rho*daz*cos(el)*cos(az);
DrhoZ=drho*sin(el)+rho*del*cos(el);

%Determine final Rho_sez and Drho_sez vectors
Rho_sez=[RhoS,RhoE,RhoZ];
Drho_sez=[DrhoS,DrhoE,DrhoZ];

%End function
end

