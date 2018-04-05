%% UpdateOrb
%%
function [ n, ecc, raan, argp, nu ] = UpdateOrb( deltat, n0, ndot2, ecc0, eccdot, raan0, raandot, argp0, argpdot, mean0 )
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%-
%
%  UpdateOrb uses the method of general perturbations to update the COEs
%  from time t0 to time t0+deltat for inclined elliptical orbits. This
%  algorithm takes into account the effects due to first-order secular
%  rates (second order for mean) caused by drag and J2 effects.
%
%  Author      : C2C Alekos Michael  USAFA  719-333-4529  15 Apr 2016
%
%  Inputs      :
%   deltat          Elapsed Time Since Epoch Date                   sec	
%   n0              Mean Motion at t0                               rad/sec
%   ndot2           Mean Motion Rate Divided by 2                 rad/sec^2
%   ecc0            Eccentricity at t0                             unitless
%   eccdot          Eccentricity Rate                               1/sec
%   raan0           RAAN at t0                                      rad
%   raandot         RAAN  rate                                      rad/sec
%   argp0           Argument of Periapsis at t0                     rad
%   argpdot         Argument of Periapsis Rate                      rad/sec
%   mean0           Mean Anomaly at t0                              rad
%
%  Outputs     :
%   n               Mean Motion at t=deltat                         rad/sec
%   ecc             Eccentricity at t = deltat                     unitless
%   raan            RAAN at t = deltat                              rad
%   argp            Argument of Periapsis at t = deltat             rad
%   nu              True Anomaly at t = deltat                      rad
%
%  Locals      :
%   E               Eccentric Anomaly at t=deltat                   rad
%   ndot0           Mean Motion Rate Multiplied by 2              rad/sec^2
%   mean            Mean Anomaly at t=deltat                        rad
%
%  Globals     :
%   None.
%
%  Coupling    :
%   newton          Solves Kepler's equation for new E
%   revcheck        Accomplishes a modulus function of x by modby.
%
%  References  :
%   None.
%
%  Documenation:
%   None.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%-

ndot0=ndot2*2; %converts TLE ndot2 to ndot0
n=n0+ndot0*deltat; %updates n
ecc=ecc0+eccdot*deltat; %updates ecc
raan=raan0+raandot*deltat; %updates raan
argp=argp0+argpdot*deltat; %updates argp
mean=mean0+n0*deltat+ndot2*deltat^2; %updates mean
mean=revcheck(mean,2*pi); %ensures mean is less than 2*pi
[E] = newton( mean,ecc ); %Calls newton to get updated eccentric anomaly 
%value
nu=acos((cos(E)-ecc)/(1-ecc*cos(E))); %updates true anomaly
if mean > pi %checks if mean is greater than 180 deg so that nu will match
    nu=2*pi-nu;
end

end

