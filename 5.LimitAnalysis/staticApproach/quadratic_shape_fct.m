function [N,DN] = quadratic_shape_fct(xi,eta)
%QUADRATIC_SHAPE_FCT Summary of this function goes here
%   Detailed explanation goes here

N(1) = xi*(2*xi-1);
N(2) = eta*(2*eta-1);
N(3) = (1-xi-eta)*(1-2*xi-2*eta);
N(4) = 4*xi*eta;
N(5) = 4*eta*(1-xi-eta);
N(6) = 4*xi*(1-xi-eta);

% compute derivatives wrt to xi
DN(1,1) = 4*xi-1;
DN(2,1) = 0;
DN(3,1) = -3+4*xi+4*eta;
DN(4,1) = 4*eta;
DN(5,1) = -4*eta;
DN(6,1) = 4-8*xi-4*eta;

% compute derivatives wrt to eta
DN(1,2) = 0;
DN(2,2) = 4*eta-1;
DN(3,2) = -3+4*xi+4*eta;
DN(4,2) = 4*xi;
DN(5,2) = 4-4*xi-8*eta;
DN(6,2) = -4*xi;
end

