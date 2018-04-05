%% PREDICT
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%-
%
%  PREDICT takes a TLE and computes the visibility conditions for the
%  satellite specified by the TLE from a given Earth location. This program
%  accomplishes this task by outputting the topocentric range, azimuth, and
%  elevation each time the Earth location is in darkness and the satellite
%  is above the location's local horizon at a reasonable range. This
%  information is calculated and output every 2 minutes (output only for
%  good visibility conditions) from a given start Julian Day to an end
%  Julian Day.
%
%  Author      : C2C Alekos Michael  USAFA  719-333-4529  15 Apr 2016
%
%  Inputs      :
%   Predict.dat     Two-Line-Element Information
%
%  Outputs     :
%   Month           Month Observation Will Occur
%   Day             Day Observation Will Occur
%   Hr              Hour Observation Will Occur (Universal Time)
%   Min             Minute Observation Will Occur (UT)
%   Rho             Range Satellite Will be Observed at             km
%   Az              Satellite Azimuth During Observation            deg
%   El              Satellite Elevation During Observation          deg
%   Vis?            Yes or No Depending on Satellite Visibility
%   a               Semi-major Axis                                 km
%   e               Eccentricity                                 unitless
%   i               Inclination                                     deg
%   raan            Right Ascension of Ascending Node               deg
%   argp            Argument of Perigee                             deg
%   nu              True Anomaly                                    deg
%   Julian Date     Days Elapsed Since 1 Jan 4713 BC
%   gst             Greenwich Sidereal Time                 0.0 to 2Pi rad
%   lst             Local Sidereal Time                     0.0 to 2Pi rad
%   deltat          Time Since Epoch                                sec
%   R               Position Vector of Satellite IJK Frame          km
%   Rsite           Position Vector of Observation Site (IJK)       km
%   RhoIJK          Range Vector of Satellite (IJK)                 km
%   RhoSEZ          Range Vector of Satellite (SEZ)                 km
%   Rsun            Position Vector of Sun (IJK)                    km
%   counter         Index of Number of Main Loop Iterations
%   Beta            Angle Between Rsun and R                        rad
%   Alpha           Angle Between Rsun and Rsite                    rad
%
%
%  Locals      :
%    fin            Open Data File for TLE Information
%    fout           Output Data File for Results
%    t0             Start Time at Start Julian Day                  sec
%    jd             Current Julian Day of Iteration
%    endtime        Time at Stop Julian Day                         sec
%    deltat         Time Since Epoch                                sec
%    i              Index of Number of Main Loop Iterations
%
%  Globals     :
%    Rad            Conversion from degrees to radians
%    MU             Earth Gravitational Constant                   km^3/s^2
%    OmegaEarth     Angular Speed of Earth                          rad/s
%
%  Coupling    :
%    input_predict  Reads TLE and converts information to working variables
%                   and units. Also outputs the epoch Julian Day
%    userinput      Prompts user for dates of observations (start and stop
%                   time as well as year). Also initializes site
%                   information
%    wgs84data      Initializes all Necessary Global Constants
%    gstlst         Computes gst and lst at respective Julian Date
%    site           Computes IJK position vector of radar site on Earth
%    rvector        Computes spacecraft IJK position vector from COEs
%    UpdateOrb      Updates the COEs given initial conditions of orbit,
%                   deltat and rates of change of COEs
%    output_predict Outputs observation information about satellites under
%                   certain visibility criterion
%    visible        Determines if current observations time meets
%                   visibility criterion.
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
wgs84data 
global Rad MU OmegaEarth

%Predict Input and Echo Checking
%Open Read File to read input data from
fin=fopen('brightest.dat','rt');
%Open Write File to Echo Check to
fout=fopen('PredictEchoCheck.out','wt');
%Define Ground Station Observation Information
[startjd,stopjd,sitlat,sitlon,sitalt] = userinput(fout);
%Create Loop to Updated COEs, calculate pertinent data/information, and 
%determine if satellite is visible.
while not (feof(fin)) %Continue loop while file not at the end of the input
    %file
    %Read input data from input_predict function and perform echo checks
    [epoch_jd,ndot2,inc0,raan0,ecc0,argp0,mean0,n0,raandot,argpdot,eccdot]...
        =input_predict(fin,fout);
    %Get epoch_jd, calculate t0 and endtime, initialize deltat, counter (i)
    %and jd
    t0=(startjd-epoch_jd)*86400;
    endtime=(stopjd-epoch_jd)*86400;
    deltat=t0;
    i=1;
    jd=startjd;

    while deltat<endtime %continues loop until stopjd reached
    %UpdateOrb updates the COEs
    [n, ecc, raan, argp, nu] = UpdateOrb( deltat, n0, ndot2, ecc0, eccdot, raan0, raandot, argp0, argpdot, mean0 );
    [R_ijk] = Rvector( n,ecc,inc0,raan,argp,nu ); %calculates Rijk
    [ gst,lst ] = gstlst( jd,sitlon ); %gets gst and lst
    R_site = site( sitlat,sitalt,lst ); %gets Rsite vector
    [rho,az,el,vis] = Visible( R_ijk,R_site,sitlat,lst,jd ); %gets Visible
    %conditions
    %loop below only outputs information if satellite is visible
   if vis == 1
       vis = 'Yes';
        output_predict(vis,rho,az,el,jd,fout);
   end
    deltat=t0+120*i; %updates deltat for next iteration
    jd=epoch_jd+deltat/86400; %updates jd for next iteration
    i=i+1; %updates loop counter
    end
%End loop

end
%Close Opened Files
fclose(fin);
fclose(fout);
