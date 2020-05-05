function [a,w] = quadrature_1D( ngauss,method)
%QUADRATURE_1D Summary of this function goes here
%   Detailed explanation goes here

switch method
    case 'gauss'
        switch ngauss
            case 1
                a = 0;
                w = 2;
            case 2
                a = [-1/sqrt(3) 1/sqrt(3)];
                w = [1 1];
            case 3   
                a = [-sqrt(3/5) 0 sqrt(3/5)];
                w = [5/9 8/9 5/9];
            case 5
                a = [-0.90617985 -0.53846931 0 0.53846931 0.90617985];
                w = [0.23692689 0.47862867 0.56888889 0.47862867 0.23692689];
            otherwise
                error('Number of edge gauss points incorrect');
        end
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

