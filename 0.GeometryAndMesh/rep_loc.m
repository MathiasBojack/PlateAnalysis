function [Rloc] = rep_loc(T)
%REP_LOC Summary of this function goes here
%   Detailed explanation goes here

T = T';
a1 = T(:,2) - T(:,1);
a2 = T(:,3) - T(:,1);


Rz = cross(a1,a2);
Rz = Rz./norm(Rz);

if (abs(Rz(2))==1)       % local z is equal to global Y
    Rx = [1;0;0];
else
    Rx = [Rz(3);0;-Rz(1)]./sqrt(Rz(1)^2+Rz(3)^2);
end

% Rx = a1./norm(a1);

Ry = cross(Rz,Rx);


Rloc = [Rx Ry Rz];
end

