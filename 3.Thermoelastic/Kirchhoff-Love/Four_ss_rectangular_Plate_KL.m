function Solution = Four_ss_rectangular_Plate_KL( Input )
%==========================================================================
%
% This function calculates the transverse displacement of a rectangular
% plate under thermal loading by using the Kirchhoff-Love plate model
% 
%=========================================================================
a   = Input.Geometry.a; % hauteur
b   = Input.Geometry.b; % largeur
CHI_SS = Input.Properties.CHI_SS; % thermal curvature

% Input.SolParam.M=50;
% Input.SolParam.N=50;
Mmax = Input.SolParam.M;
Nmax = Input.SolParam.N;
N=100;

dx=a/N;
dy=b/N;

x=0:dx:a;
y=0:dy:b;
W=zeros(N+1);
Chi_11=zeros(N+1);
Chi_22=zeros(N+1);
Chi_12=zeros(N+1);

[xx,yy]=meshgrid(x,y);
for m=1:2:Mmax
    for n=1:2:Nmax
        alpha = m*pi/a;
        beta  = n*pi/b;
        Wmn=-16*CHI_SS/pi^2/m/n/(alpha^2+beta^2);
        W=W+Wmn*sin(alpha*xx).*sin(beta*yy);
        Chi_11 = Chi_11-Wmn*alpha^2*sin(alpha*xx).*sin(beta*yy);
        Chi_22 = Chi_22-Wmn*beta^2*sin(alpha*xx).*sin(beta*yy);
        Chi_12 = Chi_12+Wmn*alpha*beta*cos(alpha*xx).*cos(beta*yy);
    end
end

Nu = Input.Properties.Nu;

Solution.Chi11=Chi_11;
Solution.Chi22=Chi_22;
Solution.Chi12=Chi_12;

Solution.W =W;
Solution.X = xx;
Solution.Y = yy; 

% surf(xx,yy,W)
end

