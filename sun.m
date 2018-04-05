function [RSun, RtAsc, Decl] = sun ( JD)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%-
%%
%%  Use           : [RSun, RtAsc, Decl] = sun ( JD)
%%
%%  This function calculates the Geocentric Equatorial position vector for
%%    the Sun given the Julian Date.  This is the low precision formula and
%%    is valid for years from 1950 to 2050.  Accuaracy of apparent coordinates
%%    is 0.01 degrees.  Notice many of the calculations are performed in
%%    degrees, and are not changed until later.  This is due to the fact that
%%    the Almanac uses degrees exclusively in their formulations.
%%
%%  Algorithm     : Calculate the several values needed to find the vector
%%                  Be careful of quadrant checks
%%
%%  Author        : Capt Dave Vallado  USAFA/DFAS  719-472-4109  25 Aug 1988
%%  In Ada        : Dr Ron Lisowski    USAFA/DFAS  719-472-4110  17 May 1995
%%  In MatLab     : Dr Ron Lisowski    USAFA/DFAS  719-333-4109  10 Oct 2001
%%
%%  Inputs        :
%%    JD          - Julian Date                            days from 4713 B.C.
%%
%%  Outputs       :
%%    RSun        - IJK Position vector of the Sun         km
%%    RtAsc       - Right Ascension                        rad
%%    Decl        - Declination                            rad
%%
%%  Locals        :
%%    MeanLong    - Mean Longitude
%%    MeanAnomaly - Mean anomaly
%%    N           - Number of days from 1 Jan 2000
%%    EclpLong    - Ecliptic longitude
%%    Obliquity   - Mean Obliquity of the Ecliptic
%%
%%  Constants     :
%%    Pi          -
%%    TwoPI       -
%%    InvRad      - Radians per degree
%%
%%  Coupling      :
%%    dfasmath.mod    Overloaded Long_Float MOD function
%%    dfasmath.ArcSin Arc Sine function
%%
%%  References             :
%%    1996 Astronomical Almanac Pg. C24
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%-

%%  Constants
    global  TwoPI Rad

     %%%%%%%%%%%%%%%%%%%%  Initialize values   %%%%%%%%%%%%%%%%%%-
     N= ( JD - 2451545.0 );

     MeanLong= 280.461 + 0.9856474*N;
     MeanLong=  revcheck (MeanLong, 360.0); %%deg%%

     MeanAnomaly= 357.528 + 0.9856003*N;
     MeanAnomaly= revcheck (MeanAnomaly*Rad, TwoPI); %%rad%%
     if MeanAnomaly < 0.0 ,
         MeanAnomaly= TwoPI + MeanAnomaly;
     end

     EclpLong= MeanLong + 1.915*sin(MeanAnomaly) + 0.020*sin(2.0*MeanAnomaly);%%deg%%
     Obliquity= 23.439 - 0.0000004*N; %%deg%%

     MeanLong = MeanLong*Rad;
     if MeanLong < 0.0 ,
         MeanLong= TwoPI + MeanLong;
     end

     EclpLong = EclpLong *Rad;
     Obliquity= Obliquity *Rad;

     RtAsc= atan( cos(Obliquity)*tan(EclpLong) );

     %%%%%% Check that RtAsc is in the same quadrant as EclpLong %%%%-

     if EclpLong < 0.0 ,
         EclpLong= EclpLong + TwoPI;   %% make sure it's in 0 to 2pi range %%
     end

     if abs( EclpLong-RtAsc ) > pi*0.5 ,
         RtAsc= RtAsc + 0.5*pi* round( (EclpLong-RtAsc)/(0.5*pi) );
     end

     Decl = asin( sin(Obliquity)*sin(EclpLong) );

     %%%%%%%%%%- Find magnitude of SUN vector, , components %%%%%%-
     RSun4= 1.00014 - 0.01671*cos( MeanAnomaly ) ...
                       - 0.00014*cos( 2.0*MeanAnomaly );   %% in AU's %%
     RSun = [RSun4*cos( EclpLong); ...
             RSun4*cos(Obliquity)*sin(EclpLong); ...
             RSun4*sin(Obliquity)*sin(EclpLong)];

     %%- Change to km from AU %%

     RSun = RSun * 149597870.0;  %% km %%
