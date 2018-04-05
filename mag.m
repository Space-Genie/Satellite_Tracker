function B=mag(A)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%                                                                     
%%  Use           : B=mag(A)
%%                                                                     
%%   Author       : Dr. RON LISOWSKI, DFAS,      5 Jan 95
%%   In MatLab    : Thomas L. Yoder, LtC, USAFA, Spring 00
%%                                                                     
%% Overview: This function Calculates the Magnitude of a Vector
%%           Using this function in lieu of the general norm function
%%            may speed execution time in large data processing tasks
%%
%%
%%   Input        :
%%      A         - Input Vector   
%%
%%   Output       :
%%      B         - Output magnitude float
%%                                                                     
%%                                                                     
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[M,N]=size(A);
if M>N  % for column vector
   B=sqrt(A(1,1)^2+A(2,1)^2+A(3,1)^2);
else    % for row vector
   B=sqrt(A(1,1)^2+A(1,2)^2+A(1,3)^2);
end



   

