%% GSTLST
%%
function [ gst,lst ] = gstlst( jd,sitlon )
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%-
%
%  Use         : [ gst,lst ] = gstlst( jd,sitlon )
%
%  The gstlst function computes the Greenwich Sidereal Time and Local 
%  Sidereal time at the time contained in the input Julian date.
%
%  Author      : C2C Alekos Michael  USAFA  719-333-4529  16 Feb 2016
%
%  Inputs      :
%    sitlon      - Site Longitude                              (rad)
%    jd          - Julian Date                                
%
%  Outputs     :
%    gst         - Greenwich Sidereal Time                0.0 to 2Pi rad
%    lst         - Local Sidereal Time                    0.0 to 2Pi rad
%
%  Locals      :
%    Dc          - Elapsed Days         
%    GST0        - gst on 1 Jan at 0000 for corresponding year
%    g           - Angle of gst since GST0   
%    l           - lst before revcheck
%
%  Globals     :
%    Rad         - Conversion from degrees to radians
%    MU          - Earth Gravitational Constant           (km^3/s^2)
%    OmegaEarth  - Angular Speed of Earth                    (rad/s)
%
%  Coupling    :
%    invjulianday- Finds year, month, day, hour, minute, and second given 
%                  the Julian date
%    finddays    - Finds fractional days through a year given the year,
%                  month, day, hour, minute and second
%    gstim0      - Finds Greenwich Sidereal time at beginning of a year.
%                  Only good for 0 hr UT, 1 Jan of a year.
%    revcheck    - Accomplishes a modulus function of x by modby.
%
%  References  :
%    None.
%
%  Documenation:
%    None.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%-
%Initialize Program Setup and Define Globals/Constants
wgs84data;
global MU
global OmegaEarth
global Rad

%Convert jd to working date and time
[Yr,Mon,D,H,M,S]=invjulianday(jd);
%Find elapsed days
Dc=finddays(Yr,Mon,D,H,M,S);
%Find gst0 on 1 Jan 0000 of calculated year
GST0=gstim0(Yr);
%Find gst since gst0
g=(GST0+1.002737791737697*(2*pi)*Dc);
%Revcheck to find true value of gst
gst=revcheck(g,2*pi);
%Calculate lst
l=gst+sitlon;
%Perform revcheck
lst=revcheck(l,2*pi);
%End function
end

