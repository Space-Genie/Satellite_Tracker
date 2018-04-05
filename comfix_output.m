function comfix_output(outfile, Rho_sez, Drho_sez, R_site, Rho_ijk, ...
                       Drho_ijk, R_ijk, V_ijk, gst, lst, A, Ecc, Incl, ...
                       RAAN, Argp, Nu)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%                       COMFIX_output
%
%   Outputs the relevant data from COMFIX
%
%   comfix_output(outfile, Rho_sez, Drho_sez, R_site, Rho_ijk, Drho_ijk,
%      R_ijk, V_ijk, gst, lst, a, ecc, incl,raan, argp, nu);
%
%   Author       : Dr Paul Vergez, USAFA, Fall 03
%
%   Inputs
%	    outfile   - File identifier for output file
%      jd        - Julian Date
%      Rho_sez   - satellite position vector from site (SEZ)         (km)
%      Drho_sez  - satellite rate vector from site (SEZ)             (km/s)
%		R_site    - vector to viewing site from Earth center (IJK)    (km)
%    	Rho_ijk   - satellite position vector from site (IJK)         (km)
%      Drho_ijk  - satellite rate vector from site (IJK)             (km/s) 
% 		R_ijk     - satellite position vector from Earth center (IJK) (km)
%      V_ijk     - satellite velocity vector from Earth center (IJK) (km/s) 
%      gst,lst   - Greenwich and local sidereal times                (rad)
%      a         - Semimajor Axis                                    (km)
%      ecc       - Eccentricity
%      incl      - Inclination                                   rad [deg]
%      raan      - Right Ascension of Ascending Node             rad [deg]
%      argp      - Argument of Perigee                           rad [deg]
%      nu        - Initial True Anomaly                          rad [deg]
%
%   Locals:
%
%   Globals:
%      Deg       - Conversion from radians to degrees
%       
%   Couplings:
%      WGS84Data - Loads global variables
%      InvJulianDay - Finds Month, Day, Hour, Min, Seconds for a given
%                     Julian day (integer or decimal) 
%      FindDays     - Converts Month, Day, Hour, Min, Seconds into a
%                     decimnal number of days since 1 Jan 00:00
%
%   References: None.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

global Deg

fprintf(outfile,'GST  = ');
fprintf(outfile,'%22.10f',gst);fprintf(outfile,' rad     &');
fprintf(outfile,'%22.10f',gst*Deg);fprintf(outfile, ' degs\n');
fprintf(outfile,'LST  = ');
fprintf(outfile,'%22.10f',lst);fprintf(outfile,' rad     &');
fprintf(outfile,'%22.10f',lst*Deg);fprintf(outfile, ' degs\n\n');
   
%%%%%%%%%%%%%%%%%%%%%%%% COMFIX output statements %%%%%%%%%%%%%%%%%%%%%%
 
fprintf(outfile,'             ####       Vectors        ####\n\n');

fprintf(outfile,'Rho_sez  = ');
fprintf(outfile,'%12.4f',Rho_sez (1)); fprintf(outfile,' S');
fprintf(outfile,'%12.4f',Rho_sez (2)); fprintf(outfile,' E');
fprintf(outfile,'%12.4f',Rho_sez (3)); fprintf(outfile,' Z');
fprintf(outfile,' Mag = '); fprintf(outfile,'%11.4f',norm(Rho_sez));
fprintf(outfile,' km \n');

fprintf(outfile,'Drho_sez = ');
fprintf(outfile,'%12.4f',Drho_sez (1)); fprintf(outfile,' S');
fprintf(outfile,'%12.4f',Drho_sez (2)); fprintf(outfile,' E');
fprintf(outfile,'%12.4f',Drho_sez (3)); fprintf(outfile,' Z');
fprintf(outfile,' Mag = '); fprintf(outfile,'%11.4f',norm(Drho_sez));
fprintf(outfile,' km/s \n');

fprintf(outfile,'R_site   = ');
fprintf(outfile,'%12.4f',R_site(1)); fprintf(outfile,' I');
fprintf(outfile,'%12.4f',R_site(2)); fprintf(outfile,' J');
fprintf(outfile,'%12.4f',R_site(3)); fprintf(outfile,' K');
fprintf(outfile,' Mag = '); fprintf(outfile,'%11.4f',norm(R_site));
fprintf(outfile,' km \n');

fprintf(outfile,'Rho_ijk  = ');
fprintf(outfile,'%12.4f',Rho_ijk(1)); fprintf(outfile,' I');
fprintf(outfile,'%12.4f',Rho_ijk(2)); fprintf(outfile,' J');
fprintf(outfile,'%12.4f',Rho_ijk(3)); fprintf(outfile,' K');
fprintf(outfile,' Mag = '); fprintf(outfile,'%11.4f',norm(Rho_ijk));
fprintf(outfile,' km \n');

fprintf(outfile,'Drho_ijk = ');
fprintf(outfile,'%12.4f',Drho_ijk(1)); fprintf(outfile,' I');
fprintf(outfile,'%12.4f',Drho_ijk(2)); fprintf(outfile,' J');
fprintf(outfile,'%12.4f',Drho_ijk(3)); fprintf(outfile,' K');
fprintf(outfile,' Mag = '); fprintf(outfile,'%11.4f',norm(Drho_ijk));
fprintf(outfile,' km/s \n');

fprintf(outfile,'R_ijk    = ');
fprintf(outfile,'%12.4f',R_ijk(1)); fprintf(outfile,' I');
fprintf(outfile,'%12.4f',R_ijk(2)); fprintf(outfile,' J');
fprintf(outfile,'%12.4f',R_ijk(3)); fprintf(outfile,' K');
fprintf(outfile,' Mag = '); fprintf(outfile,'%11.4f',norm(R_ijk));
fprintf(outfile,' km \n');

fprintf(outfile,'V_ijk    = ');
fprintf(outfile,'%12.4f',V_ijk(1)); fprintf(outfile,' I');
fprintf(outfile,'%12.4f',V_ijk(2)); fprintf(outfile,' J');
fprintf(outfile,'%12.4f',V_ijk(3)); fprintf(outfile,' K');
fprintf(outfile,' Mag = '); fprintf(outfile,'%11.4f',norm(V_ijk));
fprintf(outfile,' km/s \n');

fprintf(outfile,'------------------------------------------------------------------\n\n');
fprintf(outfile,'             ####     Output Data      ####\n\n');

fprintf(outfile,'CLASSIC ORBITAL ELEMENTS: \n\n');

fprintf(outfile,'Semimajor Axis                     (km) = ');
fprintf(outfile,'%14.4f\n',A);

fprintf(outfile,'Eccentricity                            = ');
fprintf(outfile,'%14.4f\n',Ecc);

fprintf(outfile,'Inclination                       (deg) = ');
fprintf(outfile,'%14.4f\n',Incl*Deg);

fprintf(outfile,'Right Ascension of Ascending Node (deg) = ');
fprintf(outfile,'%14.4f\n',RAAN*Deg);

fprintf(outfile,'Argument of Perigee               (deg) = ');
fprintf(outfile,'%14.4f\n',Argp*Deg);

fprintf(outfile,'True Anomaly                      (deg) = ');
fprintf(outfile,'%14.4f\n\n\n\n\n\n\n',Nu*Deg);

