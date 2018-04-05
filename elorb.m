function [P,A,Ecc,Incl,RAAN,Argp,Nu,M,U,L,CapPi] = elorb ( R,V )
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%-
%%
%%  Use           : [P,A,Ecc,Incl,RAAN,Argp,Nu,M,U,L,CapPi] = elorb ( R,V )
%%
%%  This function finds the classical orbital elements given the Geocentric
%%    Equatorial Position and Velocity vectors.  Special cases for equatorial
%%    and circular orbits are also handled.  if elements are Infinite, they
%%    are set to inf (IEEE Standard).  If elements are Undefined, they are 
%%    set to NaN (IEEE Standard).
%%    Be sure to check for these during output using isinf or isnan!!
%%
%%  Algorithm     : Initialize variables
%%                  if the HBar magnitude exists, continue, otherwise exit and
%%                       assign undefined values
%%                    Find vectors and values
%%                    Determine the type of orbit with if statements
%%                    Find angles depending on the orbit type
%%
%%  Author        : Capt Dave Vallado  USAFA/DFAS  719-472-4109  20 Sep 1990
%%  In Ada        : Dr Ron Lisowski    USAFA/DFAS  719-472-4110  17 May 1995
%%  Updated       : Dr Ron Lisowski    USAFA/DFAS  719-472-4110  2  Jan 2001
%%  In MatLab     : Dr Ron Lisowski    USAFA/DFAS  719-333-4109  10 Oct 2001
%%
%%  Inputs        :
%%    R           - IJK Position vector                    km
%%    V           - IJK Velocity vector                    km / sec
%%
%%  Outputs       :
%%    P           - Semi-latus rectum                      km
%%    A           - semi-major axis                        km
%%    Ecc         - eccentricity
%%    Incl        - inclination                            0.0 to Pi rad
%%    RAAN        - Right Ascension of Ascending Node      0.0 to 2Pi rad
%%    Argp        - Argument of Perigee                    0.0 to 2Pi rad
%%    Nu          - True anomaly                           0.0 to 2Pi rad
%%    M           - Mean anomaly                           0.0 to 2Pi rad
%%    U           - Argument of Latitude       (CI)        0.0 to 2Pi rad
%%    L           - True Longitude             (CE)        0.0 to 2Pi rad
%%    CapPi       - Longitude of Periapsis     (EE)        0.0 to 2Pi rad
%%
%%  Locals        :
%%    Hbar        - Angular Momentum H Vector              km^2 / sec
%%    Ebar        - Eccentricity     E Vector
%%    Nbar        - Line of Nodes    N Vector
%%    c1          - V**2 - u/R
%%    RDotV       - R dot V
%%    Hk          - Hk unit vector
%%    SME         - Specfic Mechanical Energy              km^2 / sec^2
%%    i           - index
%%    E           - Eccentric Anomaly                      rad
%%    D           - Parabolic Eccentric Anomaly            rad
%%    F           - Hyperbolic Eccentric Anomaly           rad
%%    Temp        - Temporary variable
%%    TypeOrbit   - Type of orbit                          EE, EI, CE, CI
%%    SinE        - Sine of Eccentric Anomaly              rad
%%    CosE        - Cosine of Eccentric Anomaly            rad
%%
%%  Constants     :
%%    HalfPi      -
%%    Pi          -
%%    TwoPI       -
%%    inf         - Flag for an infintite element
%%    NaN         - Flag for an undefined element
%%    Small       - Tolerance for roundoff errors
%%    Zero_IE     - Tolerance for zero inclination & eccentricity
%%
%%  Coupling      :
%%    vecangle         Find the angle between two vectors
%%
%%  References    :
%%    BMW           pg. 58 - 71, 181-188
%%    Escobal       pg. 104-107
%%    Kaplan        pg. 29 - 37
%%    Vallado       pg. 215 Ed. 1
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%-

%%   CONSTANTS
     global  TwoPI Small Zero_IE MU Rad

     %%%%%%%%%%%%%%%%%%%%-  Find H N and E vectors   %%%%%%%%%%%%%%%%
     Hbar = cross (R, V);
     RDotV = dot (R,V);

     if mag (Hbar) > Small ,
     %% if Hbar(Mag) > 0 the orbit is possible: Start%%
         Nbar = [-Hbar(2); Hbar(1); 0.0];
         Ebar = cross (V,Hbar) / MU - R / mag(R);


     %%%%%%%%%%%%%%%%%%  Find a, e and semi-latus rectum  %%%%%%%%%%%%
         SME = mag (V)^2 * 0.5  -  MU / mag(R);
         if abs( SME ) > Small ,
             A = -MU / (2.0*SME);       %% All Conics but Parabola %%
         else
             A= inf;                           %% Parabola %%
         end

         Ecc= mag (Ebar);
         P = mag (Hbar)^2 / MU;


     %%%%%%%%%%%%%%%%%%%%%%  Find inclination   %%%%%%%%%%%%%%%%%%%%%%
         Hk= Hbar(3)/mag(Hbar);

         if abs( abs(Hk) - 1.0 ) < Zero_IE ,
             if abs(Hbar(3)) > 0.0 ,
     %%%%%%%%%%%%           Equatorial Orbits           %%%%%%%%%%%%%%
                 Hk = sign(Hbar(3));
             end
         end
         Incl= acos( Hk );

     %%%%%%%%%%%%%%%%%%%%%%%  Find Orbit Type %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

     TypeOrbit = 'EI';
     if Incl < Zero_IE*Rad || abs(pi-Incl) < Zero_IE*Rad,
         TypeOrbit(2) = 'E';
     end
     if Ecc < Zero_IE,
         TypeOrbit(1) = 'C';
     end

     %%%%%%%%%%%%%%  Calculate 6 Angles Based on OrbitType  %%%%%%%%%%%%%%%%

     %%%%%%%%%%%%%%     Initialize All to be Undefined     %%%%%%%%%%%%%%%%
         RAAN = NaN;
         Argp = NaN;
         Nu = NaN;
         U = NaN;
         CapPi = NaN;
         L = NaN;
     
     switch upper(TypeOrbit)
     case 'EI',   % Elliptical, Parabolic or Hyperbolic Inclined

     %%%%%%%%%%%%%%  Find Right Ascension of Ascending Node %%%%%%%%%%%%%%%%
         RAAN = vecangle( Nbar, [1;0;0] );
         if isfinite(RAAN) && (Nbar(2) < 0.0) ,
            RAAN = TwoPI - RAAN;
         end
     %%%%%%%%%%%%%%%%%% Find Argument of perigee %%%%%%%%%%%%%%%%%%
         Argp = vecangle( Nbar,Ebar);
         if isfinite(Argp) && (Ebar(3) < 0.0) ,
            Argp = TwoPI - Argp;
         end
     %%%%%%%%%%%%%%  Find True Anomaly at Epoch    %%%%%%%%%%%%%%%%
         Nu = vecangle( Ebar, R );
         if isfinite(Nu) && (RDotV < 0.0) ,
            Nu= TwoPI - Nu;
         end
     %%%%%%%%%%%%%%     Others Remain Undefined    %%%%%%%%%%%%%%%%

     case 'EE',   % Elliptical, Parabolic or Hyperbolic Equatorial
     %%%%%%%%%%%%%%  Find True Anomaly at Epoch    %%%%%%%%%%%%%%%%
         Nu = vecangle( Ebar, R );
         if isfinite(Nu) && (RDotV < 0.0) ,
            Nu= TwoPI - Nu;
         end
     %%%%%%%%%%%%%%    Find Longitude of Perigee     %%%%%%%%%%%%%%
         CapPi = vecangle (Ebar, [1;0;0]);
         if isfinite(CapPi) && Ebar(2) < 0.0 ,
            U= TwoPI - U;
         end
     %%%%%%%%%%%%%%     Others Remain Undefined    %%%%%%%%%%%%%%%%
            
     case 'CI',   % Circular Inclined
     %%%%%%%%%%%  Find Right Ascension of Ascending Node %%%%%%%%%%
         RAAN = vecangle( Nbar, [1;0;0] );
         if isfinite(RAAN) && (Nbar(2) < 0.0) ,
            RAAN = TwoPI - RAAN;
         end
     %%%%%%%%%%%%%%%%  Find Argument of Latitude   %%%%%%%%%%%%%%%%
         U = vecangle( Nbar, R);
         if isfinite(U) && (R(3) < 0.0) ,
            U= TwoPI - U;
         end
     %%%%%%%%%%%%%%     Others Remain Undefined    %%%%%%%%%%%%%%%%
             
     case 'CE',   % Circular Equatorial
     %%%%%%%%%%%%%%  Find True Longitude of Epoch  %%%%%%%%%%%%%%%%
         L = vecangle (R, [1;0;0]);
         if isfinite(L) && R(2) < 0.0,
            L= TwoPI - L;
         end
     %%%%%%%%%%%%%%     Others Remain Undefined    %%%%%%%%%%%%%%%%
             
     otherwise,
          disp ('Orbit Type Error in ElOrb');
     end

     %%%%%%%%%%%%%%% Find Mean Anomaly for all orbits %%%%%%%%%%%%%%%
     % Hyperbolic and parbolic mean anomaly are not periodic angles %
     %%%%%%%%%%%%%%%%%%%%%%%%%% Hyperbolic  %%%%%%%%%%%%%%%%%%%%%%%%%
         if (Ecc-1.0) > Zero_IE ,
             F= acosh( ( Ecc+cos(Nu) ) / ( 1.0+Ecc*cos(Nu) ) );
             M= Ecc*sinh( F ) - F;
             if Nu > pi ,
                M = -M;
             end
     %%  Negative M indicates inbound, positive indicates outbound                 
         else

     %%%%%%%%%%%%%%%%%%%%%%%%%% Parabolic %%%%%%%%%%%%%%%%%%%%%%%%%%
             if abs( Ecc-1.0 ) < Zero_IE ,
                 D= sqrt( P ) * tan( Nu*0.5 );
                 M= ( 3.0*P*D + D*D*D )/6.0;
                 if Nu > pi ,
                    M = -M;
                 end
     %%  Negative M indicates inbound, positive indicates outbound                 
             else

     %%%%%%%%%%%%%%%%%%%%%%%%%% Elliptical %%%%%%%%%%%%%%%%%%%%%%%%%%
                 if Ecc > Zero_IE ,
                     Temp= 1.0 + Ecc*cos(Nu);
                     if abs(Temp) < Small ,
                         M= 0.0;
                     else
                         E = 2 * atan ( sqrt((1-Ecc)/(1+Ecc))*tan(Nu/2));
                         M = E - Ecc * sin (E);
                         if M < 0, M = TwoPI + M; end
                     end
                 else

     %%%%%%%%%%%%%%%%%%%%%%%%%%%% Circular %%%%%%%%%%%%%%%%%%%%%%%%%%
                     if  Incl < Zero_IE*Rad || abs(Incl-pi) < Zero_IE*Rad ,
                           M= L;
                     else
                           M= U;
                     end
                 end
             end
         end

      %% if Hbar(Mag) > 0 the orbit is possible: Finish%%
      else
      %% if Hbar(Mag) <= 0 the orbit is not possible: Start%%
         P    = NaN;
         A    = NaN;
         Ecc  = NaN;
         Incl = NaN;
         RAAN = NaN;
         Argp = NaN;
         Nu   = NaN;
         M    = NaN;
         U    = NaN;
         L    = NaN;
         CapPi= NaN;
      %% if Hbar(Mag) <= 0 the orbit is not possible: Finish%%
      end


