function [a,w] = quadrature_1D( ngauss,method)
%QUADRATURE_1D Summary of this function goes here
%   Detailed explanation goes here

switch method
    case 'gauss'
        [a,w] = lgwt(ngauss,-1,1);
        a = a';
        w = w';
    case 'trapeze'
        w = 2*ones(1,ngauss);
        w(1) = 1;
        w(end) = 1;
        w = w./(ngauss-1);
        a = linspace(-1,1,ngauss);
    case 'uniform'
        a = linspace(-1,1,ngauss+1);
        w = diff(a);
        a = (a(1:end-1)+a(2:end))./2;
        
    otherwise
        error('TODO');
end

end

