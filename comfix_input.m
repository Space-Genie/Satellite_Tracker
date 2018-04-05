%% COMFIX_INPUT
%%
function [sitlat,sitlon,sitalt,rho,az,el,drho,daz,del,jd] = ...
          comfix_input(fid_input, fid_output)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% COMFIX_INPUT handles inputs for COMFIX, reading/echoing input data, 
%        and converting to/echoing right units
%
% [sitlat,sitlon,sitalt,rho,az,el,drho,daz,del,jd] = ...
%          comfix_input(fid_input, fid_output)
%
% Author: Dr Paul Vergez, DFAS, Fall 2003
%         C2C Alekos Michael  USAFA  719-333-4529  16 Feb 2016
%
% Inputs:
%    fid_input  - File Identifier for input file
%    fid_output - File Identifier for output file
%
% Outputs:
%    sitlat     - Site Latitude                               (rad)
%    sitlon     - Site Longitude                              (rad)
%    sitalt     - Site Altitude                                (km)
%    rho        - Range                                        (km)
%    az         - Azimuth                                     (rad)
%    el         - Elevation                                   (rad)
%    drho       - Range Rate                                 (km/s)
%    daz        - Azimuth Rate                              (rad/s)
%    del        - Elevation Rate                            (rad/s)
%    jd         - Julian Date                                (days)
%
% Locals:
%    satnum     - Satellite Identification Number (integer) 
%    daynum     - Day of the year for observation (integer) 
%    year       - Current year (integer)
%    mon        - Month of the year for observation (integer) 
%    day        - Day of the month for observation (integer) 
%    hour       - Hour of the day for observation (integer)
%    min        - Minute of the day for observation (integer)
%    sec        - Second of the day for observation (decimal)
%    ElapsedDays- Number of days elapsed since 1 Jan, 0000 hrs (decimal)
%
% Globals:
%    Rad        - Conversion from degrees to radians
%
% Coupling:
% DayOfYr2MDHMS - Finds Month, Day, Hours, Min, and Seconds given a
%                   decimal number of days in a year
%    JulianDay  - Finds Julian Day given DayOfYr2MDHMS info
%    Finddays   - Finds fractional days through a year given the year,
%                   month, day, hour, minute and second
%       
% References   : None.
%
% Documentation: None.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

global Rad

% Read in first line of data with all the site location and time 
% information:

sitlat=fscanf(fid_input,'%f',1);
sitlon=fscanf(fid_input,'%f',1);
sitalt=fscanf(fid_input,'%f',1);
year=2000+str2double(fscanf(fid_input,'%3c',1));
daynum=fscanf(fid_input,'%3d',1);
hour  =fscanf(fid_input,'%2d',1);
min   =fscanf(fid_input,'%2d',1);
sec   =fscanf(fid_input,'%5f',1);

% Read in position and rate data

satnum=fscanf(fid_input,'%4i',1);
rho =fscanf(fid_input,'%f',1);
az  =fscanf(fid_input,'%f',1);
el  =fscanf(fid_input,'%f',1);
drho=fscanf(fid_input,'%f',1);
daz =fscanf(fid_input,'%f',1);
del =fscanf(fid_input,'%f',1);
   
% Echo-check the data

fprintf(fid_output,'*********************Comfix  Satellite: ');
fprintf(fid_output,'%4i*******************\n',satnum);
fprintf(fid_output,'---------------------------------------------------------------\n');
fprintf(fid_output,'             ####     Input Data      ####\n\n');
fprintf(fid_output,'   LAT      LON     ALT      YEAR     DAY          UT     \n');
fprintf(fid_output,'  (deg)    (deg)    (m)              Number   (hr:min:sec)\n');
fprintf(fid_output,'%7.2f',sitlat);
fprintf(fid_output,'%9.2f',sitlon);
fprintf(fid_output,'%8.2f',sitalt);
fprintf(fid_output,'%9i',year);
fprintf(fid_output,'%7i',daynum);
fprintf(fid_output,'%9i',hour);
fprintf(fid_output,':');
fprintf(fid_output,'%2i',min);fprintf(fid_output,':');
fprintf(fid_output,'%4.2f\n\n',sec);
fprintf(fid_output,'   RHO         DRHO        EL           DEL         AZ          DAZ  \n');
fprintf(fid_output,'   (km)       (km/s)     (deg)        (deg/s)     (deg)       (deg/s) \n');
fprintf(fid_output,'%9.4f',rho);
fprintf(fid_output,'%12.5f',drho);
fprintf(fid_output,'%12.5f',el);
fprintf(fid_output,'%12.5f',del);
fprintf(fid_output,'%12.5f',az);
fprintf(fid_output,'%12.5f\n',daz);

% Find out number of elapsed days, specific day of the year, and Julian date

%    >>>>>>>>>>>>> YOU FILL IN HERE! <<<<<<<<<<<<<<<<<<<
%Julian Date
[Mon,D]=dayofyr2mdhms(year,daynum); %Specific Day Of the Year
jd=julianday(year,Mon,D,hour,min,sec);
%Elapsed Days
ElapsedDays=finddays(year,Mon,D,hour,min,sec);

% Converting to preferred units

%    >>>>>>>>>>>>> YOU FILL IN HERE! <<<<<<<<<<<<<<<<<<<
sitlat=sitlat*Rad;
sitlon=sitlon*Rad;
sitalt=sitalt/1000;
el=el*Rad;
del=del*Rad;
az=az*Rad;
daz=daz*Rad;

%%%%%%%%%%%%%%%%%%%%%%%% Echo Check new units %%%%%%%%%%%%%%%%%%%%%%

fprintf(fid_output,'------------------------------------------------------------\n');
fprintf(fid_output,'             ####     Working Data    ####\n\n');
fprintf(fid_output,'     LAT        LON         ALT         Julian Date \n');
fprintf(fid_output,'    (rad)      (rad)       (km)           (days)     \n');
fprintf(fid_output,'%9.4f',sitlat);
fprintf(fid_output,'%12.4f',sitlon);
fprintf(fid_output,'%12.4f',sitalt);
fprintf(fid_output,'%21.7f\n\n',jd);
   
fprintf(fid_output,'   RHO         DRHO        EL           DEL         AZ          DAZ  \n');
fprintf(fid_output,'   (km)       (km/s)     (rad)        (rad/s)     (rad)       (rad/s) \n');
fprintf(fid_output,'%9.4f',rho);
fprintf(fid_output,'%12.5f',drho);
fprintf(fid_output,'%10.5f',el);
fprintf(fid_output,'%14.5f',del);
fprintf(fid_output,'%12.5f',az);
fprintf(fid_output,'%12.5f\n\n',daz);
