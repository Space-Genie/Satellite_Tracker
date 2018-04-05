%% COMFIX
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%-
%
%  COMFIX converts radar observation data into Classical Orbital 
%  Elements. This file also outputs the position and velocity vectors 
%  associated with the spacecraft as well as other useful spacecraft and
%  radar site information.
%
%  Author      : C2C Alekos Michael  USAFA  719-333-4529  16 Feb 2016
%
%  Inputs      :
%    sitlat      - Site Latitude                               (deg)
%    sitlon      - Site Longitude                              (deg)
%    sitalt      - Site Altitude                                 (m)
%    year        - Four Digit Year                             
%    daynum      - Julian Day of Year    
%    UT          - Observation Date in Universal Time      (HH:MM:SS.SS) 
%    satnum      - Four Digist Satellite ID Number
%    az          - Azimuth                                     (deg)
%    el          - Elevation                                   (deg)
%    rho         - Range                                        (km)
%    drho        - Range Rate                                 (km/s)
%    daz         - Azimuth Rate                              (deg/s)
%    del         - Elevation Rate                            (deg/s)
%
%  Outputs     :
%    A           - semi-major axis                              (km)
%    Ecc         - eccentricity
%    Incl        - inclination                            0.0 to 180 deg
%    RAAN        - Right Ascension of Ascending Node      0.0 to 360 deg
%    Argp        - Argument of Perigee                    0.0 to 360 deg
%    Nu          - True anomaly                           0.0 to 360 deg
%    sitlat      - Site Latitude                               (rad)
%    sitlon      - Site Longitude                              (rad)
%    sitalt      - Site Altitude                                (km)
%    rho         - Range                                        (km)
%    az          - Azimuth                                     (rad)
%    el          - Elevation                                   (rad)
%    drho        - Range Rate                                 (km/s)
%    daz         - Azimuth Rate                              (rad/s)
%    del         - Elevation Rate                            (rad/s)
%    jd          - Julian Date                               
%    Rho_sez     - Relative Position Vector in SEZ frame        (km)
%    Drho_sez    - Relative Velocity Vector in SEZ frame      (km/s)
%    R_site      - Site Earth Centered Inertial position vector (km)
%    Rho_ijk     - Relative Position Vector in IJK frame		(km)
%    Drho_ijk    - Relative Velocity Vector in IJK frame	  (km/s)
%    R_ijk       - Position Vector in IJK frame                 (km)
%    V_ijk       - Velocity vector in IJK frame               (km/s)
%    gst         - Greenwich Sidereal Time                0.0 to 2Pi rad
%    lst         - Local Sidereal Time                    0.0 to 2Pi rad
%
%  Locals      :
%    fin         - Open Data File for Radar Observation Information          
%    fout        - Output Data File for Unit Conversions and Echo Check 
%    fin2        - Reopen Data File for Radar Observation Information             
%    fout2       - Output Data File for Final Answers
%    OE          - OmegaEarth in 3D vector form in IJK frame (rad/s)
%
%  Globals     :
%    Rad         - Conversion from degrees to radians
%    MU          - Earth Gravitational Constant           (km^3/s^2)
%    OmegaEarth  - Angular Speed of Earth                    (rad/s)
%
%  Coupling    :
%  comfix_input  - Handles Inputs for COMFIX, reading/echoing input data, 
%                  and converting to and echoing correct units
%    wgs84data   - Initializes all Necessary Global Constants
%    gstlst      - Computes gst and lst at respective Julian Date
%    site        - Computes IJK position vector of radar site on Earth
%    rvtopos     - Converts radar measurement to range and range rate 
%                  vectors in SEZ frame
%    axisrot     - Performs rotation of angle ALPHA about desired axis
%    elorb       - Finds COE's given the Geocentric Equatorial Position 
%                  and Velocity vectors
%  comfix_output - Outputs the relevant data from COMFIX
%
%  References  :
%    None.
%
%  Documenation:
%    None.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%-
%Initialize Program Setup and Define Globals/Constants
clear all %clears any previously defind variables
close all %closes any open figures
clc %clears command window data
%Define Globals/Constants
wgs84data 
global Rad
global MU
global OmegaEarth
OE=[0,0,OmegaEarth].'; %turns OmegaEarth into 3D vector in IJK coordinates

%Comfix Input and Echo Checking
%Open Read File to read input data from
fin=fopen('comfix.dat','rt');
%Open Write File to Echo Check to
fout=fopen('EchoCheck.out','wt');
%Create Loop to read in observation data from multiple radar observations
while not (feof(fin)) %Continue loop while file not at the end of the file
    %Read input data from comfix_input function and perform 2 echo checks
    [sitlat,sitlon,sitalt,rho,az,el,drho,daz,del,jd] = ...
          comfix_input(fin, fout);
%End loop
end
%Close Opened Files
fclose(fin);
fclose(fout);

%Output Calculated Data and Final Answers
%Reopen Read File to Read input data from
fin2=fopen('comfix.dat','rt');
%Open Write File to Output final answers to
fout2=fopen('Answers.out','wt');
%Run Loop to read in observation data from multiple radar observations
while not (feof(fin2))
    [sitlat,sitlon,sitalt,rho,az,el,drho,daz,del,jd] = ...
          comfix_input(fin2, fout2); %read in input data again
      [gst,lst]=gstlst(jd,sitlon); %calculate gst and lst 
      R_site = site( sitlat,sitalt,lst ); %determine IJK vector of site
      R_site=R_site.'; %Turn from row vector to column vector to match 
      %matrix dimensions of other vectors in future operations
      [ Rho_sez,Drho_sez ] = rvtopos( rho,az,el,drho,daz,del ); %calculate
      %range and range rate vectors in SEZ frame
      colat=pi/2-sitlat; %Determine colatitude 
      %Rotate significant vectors from SEZ to IJK frame
    temp1=axisrot(Rho_sez,2,-colat);
    temp2=axisrot(Drho_sez,2,-colat);
    Rho_ijk=axisrot(temp1,3,-lst); 
    Drho_ijk=axisrot(temp2,3,-lst); 
    %Calculate R_ijk and V_ijk vectors
    R_ijk=R_site+Rho_ijk; 
    V_ijk=Drho_ijk+cross(OE,R_ijk); 
    %Determine the COE's using elorb
    [P,A,Ecc,Incl,RAAN,Argp,Nu] = elorb ( R_ijk,V_ijk );
    %Output final answers using comfix_output
    comfix_output(fout2, Rho_sez, Drho_sez, R_site, Rho_ijk, ...
                       Drho_ijk, R_ijk, V_ijk, gst, lst, A, Ecc, Incl, ...
                       RAAN, Argp, Nu);
%Close loop
end
%Close Opened Files
fclose(fin2);
fclose(fout2);






