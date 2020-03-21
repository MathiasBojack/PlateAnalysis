% ============itr============================================================
%                 Material properties caculated from Safir                
% ========================================================================
 function [W,K1,K2,K4,C,S,CHI]=Newton_Raphson(Input)

M       = Input.SolParam.M;
N       = Input.SolParam.N;
MN      = M*N;

a       = Input.Geometry.a;
b       = Input.Geometry.b;
h       = Input.Geometry.h;

A       = Input.Properties.A;

nu      = Input.Properties.Nu;

gamma   = Input.Load.gamma;
De      = Input.Properties.De;

c4      = gamma*a^3/De;
c_phi   = A*(1-nu^2)*h^2/De*a^4/b^4;
chi     = Input.Properties.CHI_SS;

 [K1,K2,K4,C,S,CHI] = Von_Karman(M,N,a,b,h);
tol=10^(-5);
maxval = 10000.0;
disp('Loop start')
W0      = sparse(MN,1);
for i=1:1:1000
    disp(['i = ', num2str(i)]);
    [FW,D_FW_W] =VK_Cubic_Function(W0,c4,c_phi,chi,K1,K2,K4,C,S,CHI);
    W           = W0 - D_FW_W\FW;
%     Temp        = W0;
    W0          = W;
    err=norm(FW,inf);
    if err<tol
        disp('The convergence has been achieved');
        return;
    end
    if norm(FW,inf)>maxval
        W=zeros(MN,1);
        disp('Solution diverges');
    end
end
% W=sparse(W);

