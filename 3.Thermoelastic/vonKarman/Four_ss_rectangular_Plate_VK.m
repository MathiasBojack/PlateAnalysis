function Solution = Four_ss_rectangular_Plate_VK(Input)

%==========================================================================
%
% This function calculates the transverse displacement of a rectangular
% plate under thermal loading by using the Kirchhoff-Love plate model
% 
%=========================================================================

a    = Input.Geometry.a;
b    = Input.Geometry.b;
h    = Input.Geometry.h;

% Add temperature dependent elastic properties as attributes of Input param
Input = dim_Prop_Cal( Input);

M=20;
N=20;
MN=M*N;
TEMP1=zeros(MN,1);
nu=0.2;

A       = Input.Properties.A;
NT      = Input.Properties.NT;
% ========================================================================
%   
%                           Calculer W
% 
% ========================================================================
[W,K1,K2,K4,C,S,CHI]=Newton_Raphson(Input);
% profile viewer        
% ========================================================================
%   
%                           Calculer PHI
% 
% ========================================================================

for p=1:1:MN
TEMP1(p)=W'*C{p}*W;
end
TEMP1=K1\TEMP1;
TEMP1=sparse(TEMP1) ;
PHI=A*(1-nu^2)*a*h^2/NT/b^3*TEMP1;
disp('PHI has been calculated')
% ========================================================================
%   
%                           Calculate the maximum of displacement
% 
% ========================================================================
% disp('Starts calculate the field')
x=(0:a/100:a)';
y=(0:b/100:b)';
[xx,yy]=meshgrid(x,y);
n=length(xx);
WXY=zeros(n);
% PHIXY=zeros(n);
% N11=zeros(n);
% N22=zeros(n);
% N12=zeros(n);
% % ------------------------------------------------------------------------
% %       N11 N22 N12 The membrane stresses 
% %       CHI11 CHI22 CHI12 The curvature of the surface
% %       WXY= sum_ij(wij*sin(i*pi*x/a)sin(j*pi*y/b))
% %       W=WXY+W_etoile;
% %       W_etoile=chi_th/2*(x^2-a*x+y^2-b*y)+sum_k(ak*sin(k*pi*x/a)+bk*sin(k*pi*y/b))
% %       M11 M22 M12 The bending moments 
% %       CHI11=D2W/DX2 CHI22=D2W/DY2 CHI12=D2W/DXDY
% ------------------------------------------------------------------------
% disp('Post treatment')
% disp('Step 1. Calculate boudary terms')
% CHI11_etoile=1;
% CHI22_etoile=1;
% parfor k=1:400,
%      a_2k_1=4*a^2/(2*k-1)^3/pi^3;
%      b_2k_1=4*b^2/(2*k-1)^3/pi^3;
%      CHI11_etoile=CHI11_etoile-a_2k_1*sin(pi*(2*k-1)/a*xx)*(pi*(2*k-1)/a)^2;
%      CHI22_etoile=CHI22_etoile-b_2k_1*sin(pi*(2*k-1)/b*yy)*(pi*(2*k-1)/b)^2;
% end
% CHI11_boundary=CHI11_etoile*chi_th;
% CHI22_boundary=CHI22_etoile*chi_th;
% CHI12_boundary=zeros(n);
% ------------------------------------------------------------------------
% disp('Step 2. Calculate interior field terms')
% CHI11_Interior=zeros(n);
% CHI22_Interior=zeros(n);
% CHI12_Interior=zeros(n);
for p=1:MN
    j=ceil(p/M);
    i=p-(j-1)*M;
    WXY=WXY+W(p)*h*sin(i*pi*xx/a).*sin(j*pi*yy/b);
%     PHIXY=PHIXY+NT*a*b*PHI(p)*(1-cos(2*i*pi*xx/a)).*(1-cos(2*j*pi*yy/b));
%     N11=N11+NT*a*b*PHI(p)*(1-cos(2*i*pi*xx/a)).*cos(2*j*pi*yy/b)*(2*j*pi/b)^2;
%     N22=N22+NT*a*b*PHI(p)*cos(2*i*pi*xx/a).*(1-cos(2*j*pi*yy/b))*(2*i*pi/a)^2;
%     N12=N12+NT*a*b*PHI(p)*sin(2*i*pi*xx/a).*sin(2*j*pi*yy/b)*(2*i*pi/a)*(2*j*pi/b);
%     CHI11_Interior = CHI11_Interior-W(p)*h*sin(i*pi*xx/a).*sin(j*pi*yy/b)*(i*pi/a)^2;
%     CHI22_Interior = CHI22_Interior-W(p)*h*sin(i*pi*xx/a).*sin(j*pi*yy/b)*(j*pi/b)^2;
%     CHI12_Interior = CHI12_Interior+W(p)*h*cos(i*pi*xx/a).*cos(j*pi*yy/b)*(i*pi/a)*(j*pi/b);
end
% % ------------------------------------------------------------------------
% disp('Step 3. Total field calculation ')
% CHI11=CHI11_Interior+CHI11_boundary;
% CHI22=CHI22_Interior+CHI22_boundary;
% CHI12=CHI12_Interior+CHI12_boundary;
% 
% M11_Interior=-B/A*(N11-NT)+D_eff*(CHI11_Interior+nu*CHI22_Interior)+MT;
% M22_Interior=-B/A*(N22-NT)+D_eff*(CHI22_Interior+nu*CHI11_Interior)+MT;
% 
% M11_boundary=D_eff*(CHI11_boundary+nu*CHI22_boundary);
% M22_boundary=D_eff*(CHI22_boundary+nu*CHI11_boundary);
% 
% N11=N11+gamma*(xx-a);
% M11=-B/A*(N11-NT)+D_eff*(CHI11+nu*CHI22)+MT;
% M22=-B/A*(N22-NT)+D_eff*(CHI22+nu*CHI11)+MT;
% M12=-B/A*N12+D_eff*(1-nu)*CHI12;

% profile viewer
% ========================================================================
%   
%                           Solution plot of Airy function,Membrane
%                           stresses and bending moments
% 
% ========================================================================
% 
[min_WXY,location]=min(WXY(:));
Solution.W = WXY;
Solution.X = xx;
Solution.Y = yy;
% [R_min,C_min]=ind2sub(size(WXY),location);
% 
% figure(1);surf(x,y,WXY);shading interp
% set(gcf,'color','white')
% xlabel('$x_1$','Interpreter','LaTex','FontSize',14) 
% ylabel('$x_2$','Interpreter','LaTex','FontSize',14)
% zlabel('$x_3$','Interpreter','LaTex','FontSize',14)
% 
% figure(2);surf(x,y,PHIXY/NT/a/b);shading interp
% set(gcf,'color','white')
% xlabel('$x_1$','Interpreter','LaTex','FontSize',14) 
% ylabel('$x_2$','Interpreter','LaTex','FontSize',14)
% zlabel('$x_3$','Interpreter','LaTex','FontSize',14)
% 
% figure(3);surf(x,y,N11/abs(NT));shading interp
% set(gcf,'color','white')
% xlabel('$x_1$','Interpreter','LaTex','FontSize',14) 
% ylabel('$x_2$','Interpreter','LaTex','FontSize',14)
% zlabel('$x_3$','Interpreter','LaTex','FontSize',14)
% % 
% figure(4);surf(x,y,N22/abs(NT));shading interp
% set(gcf,'color','white')
% xlabel('$x_1$','Interpreter','LaTex','FontSize',14) 
% ylabel('$x_2$','Interpreter','LaTex','FontSize',14)
% zlabel('$x_3$','Interpreter','LaTex','FontSize',14)
% % 
% figure(5);surf(x,y,N12/abs(NT));shading interp
% set(gcf,'color','white')
% xlabel('$x_1$','Interpreter','LaTex','FontSize',14) 
% ylabel('$x_2$','Interpreter','LaTex','FontSize',14)
% zlabel('$x_3$','Interpreter','LaTex','FontSize',14)
% % 
% figure(6);surf(x,y,M11/MT);shading interp
% set(gcf,'color','white')
% xlabel('$x_1$','Interpreter','LaTex','FontSize',14) 
% ylabel('$x_2$','Interpreter','LaTex','FontSize',14)
% zlabel('$x_3$','Interpreter','LaTex','FontSize',14)
% 
% figure(7);surf(x,y,M22/MT);shading interp
% set(gcf,'color','white')
% xlabel('$x_1$','Interpreter','LaTex','FontSize',14) 
% ylabel('$x_2$','Interpreter','LaTex','FontSize',14)
% zlabel('$x_3$','Interpreter','LaTex','FontSize',14)
% 
% figure(8);surf(x,y,M12/MT);shading interp
% set(gcf,'color','white')
% xlabel('$x_1$','Interpreter','LaTex','FontSize',14) 
% ylabel('$x_2$','Interpreter','LaTex','FontSize',14)
% zlabel('$x_3$','Interpreter','LaTex','FontSize',14)




        
        
