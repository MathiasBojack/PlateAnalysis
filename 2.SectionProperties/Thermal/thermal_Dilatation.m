function [ epsilon ] = thermal_Dilatation( T ,option)
%This function caculates the dilatation of concret under temperature T
%   Two different formation are considered here for the concret:
%   - Case 1:  concrete with sand
%   - Case 2:  concrete with limestone

n=length(T);
epsilon=zeros(n,1);
if nargin <2
    option = 'Silicieux';
end
switch option
    case 'Silicieux'
       Index_1 = T<=700;
       Index_2 = T>700;
       epsilon(Index_1) = -1.8e-4+9e-6*T(Index_1)+2.3e-11*T(Index_1).^3;
       epsilon(Index_2) = 14e-3;
    case 'Calcaire'
        Index_1 = T<=805;
       Index_2 = T>805;
       epsilon(Index_1) = -1.2e-4+6e-6*T(Index_1)+1.4e-11*T(Index_1).^3;
       epsilon(Index_2) = 12e-3;
end
end

