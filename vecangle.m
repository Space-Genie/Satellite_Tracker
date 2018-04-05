function Theta = vecangle(A,B)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%%  Use           : Theta = vecangle(A,B)
%%
%%  This function calculates the angle between two vectors.  The output is
%%    set to a large number to indicate an undefined value.  Be SURE to check for
%%    this at the output phase.
%%
%%  Algorithm     : Check the denominator for a divide by zero
%%                  Check for exactly 1.0 or -1.0 to avoid ArcCosine problems
%%
%%  Author        : Capt Dave Vallado     USAFA/DFAS  719-472-4109   14 Sep 1990
%%  In Ada        : Dr Ron Lisowski       USAFA/DFAS  719-472-3315   23 Jan 1995
%% In MatLab      : Thomas L. Yoder, LtC, USAFA,                     Spring 00
%%
%%  Inputs        :
%%       A         % Input Vector                    array of dimension 3
%%       B         % Input Vector                    array of dimension 3
%%
%%  OutPuts       :
%%    Theta        - Angle between the two vectors               -Pi to Pi
%%
%%  Locals        :
%%    temp        - Temporary REAL variable
%%
%%  Coupling      :
%%    dot         : Finds the dot product of two vectors
%%	  mag	      : Finds the vector magnitude	
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

Tolerance=0.000001;

if mag(A)*mag(B)>=Tolerance,
    temp=dot(A,B)/mag(A)/mag(B);
    if abs(temp)>1.0,
       temp=sign(temp)*1.0;
    end
    Theta=acos(temp);
 else
    Theta=NaN;
 end