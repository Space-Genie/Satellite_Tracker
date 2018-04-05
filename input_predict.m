%% Input_Predict
%%
function [epoch_jd,ndot2,inc0,raan0,ecc0,argp0,mean0,n0,raandot,argpdot,eccdot]=input_predict(infile,outfile)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%                       Program input_predict
%
% [epoch_jd,ndot2,inc0,raan0,ecc0,argp0,mean0,n0,raandot,argpdot,
%                            eccdot]=input_predict(infile,outfile)
%
% Given the Two-line element set for a satellite at a given epoch time, 
% this program calculates viewing opportunities and directions for the
% satellite from a given site over a time interval
%
%   Author       : Dr Paul Vergez, USAFA, Spring 2013
%
%   Input        :
%      infile       - flag for input file for TLES
%      outfile      - flag for output file
%
%   Outputs       :
%	   epochjd      - Julian Date at epoch time (when COEs known)  (days)
%      n0           - initial mean motion                         (rad/s)
%      ndot2        - mean motion divided by 2                  (rad/s^2)
%      ecc0         - initial eccentricity
%      eccdot       - rate change in eccentricity                  (/sec)
%      raan0        - initial right ascension of the asc node       (rad)
%      raandot      - rate change in RAAN                         (rad/s)
%      argp0        - initial argument of perigee                   (rad)
%      argpdot      - rate change in argp                         (rad/s)
%      mean0        - initial Mean Anomaly                          (rad)
%      inc0         - initial inclination                           (rad)
%
%   Locals       :
%		epochyear    - year at epoch time (when COEs known)
%      epochday     - day of year (decimal) at epoch time
%      satnum       - satellite identification number
%      satname      - name of satellite
%   	month        - month at epoch time (integer) 
%  	day          - Day of the month at epoch time (integer) 
%      hour         - Hour of the day at epoch time (integer)
%      min          - Minute of the hour at epoch time (integer)
%      sec          - Second of the minute at epoch time (decimal)
%
%   Constants
%      TwoPI        - 2*pi
%      Rad          - degrees to radians conversion
%
%   Coupling:
%      J2DragPert   - Calculates secular rates for ecc, raan, and argp
%                          at epoch time
%      DayOfYr2MDHMS- Finds Month, Day, Hour, Min, Seconds for a given
%							day (integer or decimal) 
%      JulianDay    - Finds Julian Day given DayOfYr2MDHMS info
%
%   References   : None.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

global TwoPI Rad

%%%%%%%%%%%%%%%%%%%%% Read data from TLE input file %%%%%%%%%%%%%%%

satname=fgetl(infile);
disp(satname)

satnum=fscanf(infile,'%*d %6s',1);
epochyear=fscanf(infile,'%*s %2d',1);

epochyear=2000+epochyear; % Change Epoch year to 4-digit #

epochday=fscanf(infile,'%f',1);

ndot2=fscanf(infile,'%f',1);
inc0=fscanf(infile,'%*s %*s %*s %*s %*s %*s %f',1);
raan0=fscanf(infile,'%f',1);
ecc0=fscanf(infile,'%d',1);ecc0=ecc0/1.0E7;
argp0=fscanf(infile,'%f',1);
mean0=fscanf(infile,'%f',1);
n0=fscanf(infile,'%f',1);
trash=fscanf(infile,'%f\n',1);

fprintf(outfile,'\n\n*********************************************************************\n');
fprintf(outfile,' Echo Check Input Data for: ');fprintf(outfile,'%15s\n',satname);
fprintf(outfile,'*********************************************************************\n\n');
fprintf(outfile,'    Epoch Year     = ');fprintf(outfile,'%12i\n',epochyear);
fprintf(outfile,'    Epoch Day      = ');fprintf(outfile,'%12.8f\n',epochday);
fprintf(outfile,'    N Dot/2        = ');fprintf(outfile,'%12.8f  rev/day^2\n',ndot2);
fprintf(outfile,'    Inclination    = ');fprintf(outfile,'%12.8f  deg\n',inc0);
fprintf(outfile,'    RAAN           = ');fprintf(outfile,'%12.8f  deg\n',raan0);
fprintf(outfile,'    Eccentricity   = ');fprintf(outfile,'%12.8f\n',ecc0);
fprintf(outfile,'    Arg of Perigee = ');fprintf(outfile,'%12.8f  deg\n',argp0);
fprintf(outfile,'    Mean Anomaly   = ');fprintf(outfile,'%12.8f  deg\n',mean0);
fprintf(outfile,'    Mean Motion    = ');fprintf(outfile,'%12.8f  rev/day\n\n',n0);

% Convert to working units:

ndot2=ndot2*(2*pi)/86400^2; %ndot2 from rev/day^2 to rad/s^2
n0=n0*(2*pi)/86400; %Mean Motion from rev/day to rad/s
inc0=inc0*Rad;
raan0=raan0*Rad;
argp0=argp0*Rad;
mean0=mean0*Rad;

% Calculate perturbations in RAAN, ArgP, and Ecc

[raandot, argpdot, eccdot] = J2DragPert(inc0, ecc0, n0, ndot2);

% Calculate Epoch Julian Day

[month,day,hour,min,sec]=dayofyr2mdhms(epochyear,epochday);
epoch_jd=julianday(epochyear,month,day,hour,min,sec);

fprintf(outfile,' Echo Check Working Units\n\n');
fprintf(outfile,'    Epoch Time (month/day, hr:min:sec) = ');
fprintf(outfile,'%3i',month);fprintf(outfile,'/');fprintf(outfile,'%2i,',day);
fprintf(outfile,'%4i',hour);fprintf(outfile,':');
fprintf(outfile,'%2i',min);fprintf(outfile,':');
fprintf(outfile,'%5.2f\n\n',sec);
fprintf(outfile,'    Epoch Julian Day (days) = ');fprintf(outfile,'%17.8f\n\n',epoch_jd);

fprintf(outfile,'    N    = %13.11f rad/s ',n0);
fprintf(outfile,'Ndot2     = %+19.11e rad/s^2\n',ndot2);
fprintf(outfile,'    Ecc  = %13.11f       ',ecc0);
fprintf(outfile,'Ecc dot   = %+19.11e 1/s \n',eccdot);
fprintf(outfile,'    Inc  = %13.11f rad   ',inc0);
fprintf(outfile,'Mean Anom = %+19.11e rad \n',mean0);
fprintf(outfile,'    RAAN = %13.11f rad   ',raan0);
fprintf(outfile,'RAAN dot  = %+19.11e rad/s\n',raandot);
fprintf(outfile,'    Argp = %13.11f rad   ',argp0);
fprintf(outfile,'Argp dot  = %+19.11e rad/s\n\n',argpdot);

%  display header for computed information to follow
fprintf(outfile,' MONTH/DAY   HR:MIN (UT)   RHO(KM)     AZ(DEG)    EL(DEG)   VIS?\n');
fprintf(outfile,' ---------------------------------------------------------------\n\n');