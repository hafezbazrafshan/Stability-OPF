function [ ZDOT ] = numericalDerivative( Z,ts )
%NUMERICALDERIVATIVE finds the numerical time derivative 
%   [ ZDOT ] = numericalDerivative( Z,ts )  calculates the numerical
%   derivative of matrix Z, where the columns show the evolution of time. 
% 





ZDOT=(Z(:,2:end)-Z(:,1:end-1))./(ts); 

end

