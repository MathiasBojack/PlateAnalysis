function [kc,ks,kt] = fire_red_coeff(Temp)
%FIRE_RED_COEFF Summary of this function goes here
%   Detailed explanation goes here


Tc = [20          50         100         200         250         300         400         500         600 ...
    700         800         900        1000        1100        1200 1e10];

kc_sampl = [1.0000    1.0000    1.0000    0.9500    0.9000    0.8500    0.7500    0.6000    0.4500    0.3000    0.1500 ...
0.0800    0.0400    0.0100         0 0];

Ts = [20         100         200         300         400         500         600         700         800 ...
    900        1000        1100        1200 1e10];

ks_sampl = [1.0000    1.0000    1.0000    1.0000    1.0000    0.7800    0.4700    0.2300    0.1100    0.0600    0.0400 ...
    0.0200         0 0];

Tt = [20         100         200         300         400         500         600         700         800         900 ...
    1000 1100 1200 1e10];
kt_sampl = [1.0000    1.0000    0.8000    0.6000    0.4000    0.2000         0         0         0         0         0         0 ...
0 0];                

kc = interp1(Tc,kc_sampl,Temp);
ks = interp1(Ts,ks_sampl,Temp);
kt = interp1(Tt,kt_sampl,Temp);
end

