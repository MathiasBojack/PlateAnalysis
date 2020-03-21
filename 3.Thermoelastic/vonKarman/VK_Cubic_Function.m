% ========================================================================
% 
%   This function caculates the cubic vetor function and the Jacobian 
%   Matrix for system of nonlinear equations, which derive from the weak 
%   form of von Karman Plate
%       
%           2016 - 10 -22 by Mingguan YANG
% =======================================================================
function [FW,D_FW_W]=VK_Cubic_Function(W,c4,c_phi,chi,K1,K2,K4,C_cell,S_cell,CHI)
% function [FW,D_FW_W]=VK_Cubic_Function(W,K1,K2,K4,C,S,CHI,gamma,a,b,A,B,D,NT,MT)
% c4    = gamma*a^3/D_eff
% c_phi = A*(1-nu^2)*h^2/D_eff*a^4/b^4

% [ A,B,D,NT,MT,~,~] = Dimensional_Sec_Cons( Time );
% [K1,K2,K4,C_array,S_array,CHI]=Von_Karman(M,N,a,b);

% h=0.15;
% nu=0.2;
% chi=-1/(A*D-B^2)*(B*NT+A*MT);
% D_eff=(A*D-B^2)/A;

MN=length(W);
TEMP1=zeros(MN,1);
TEMP2=zeros(MN,1);
D_TEMP1_DW=zeros(MN,MN);
D_TEMP2_DW=zeros(MN,MN);
% ------------------------------------------------------------------------
%                           [PHI]=A*(1-nu^2)*a*h^2/NT/b^3*[TEMP1]
%                                              {  [W]T*[C_11]*[W]  }
%                                              {         .         }
%                                              {         .         }
%                           [TEMP1]= [K1]^(-1)*{  [W]T*[C_ij]*[W]  }
%                                              {         .         }
%                                              {         .         }
%                                              {  [W]T*[C_MN]*[W]  }
%                               [y]=[TEMP1]
% ------------------------------------------------------------------------ 
%                                              {[TEMP1]T*[S_11]*[W]}
%                                              {         .         }
%                                              {         .         }
%                           [TEMP2]=           {[TEMP1]T*[S_ij]*[W]}
%                                              {         .         }
%                                              {         .         }
%                                              {[TEMP1]T*[S_MN]*[W]}
% ------------------------------------------------------------------------ 
%   F     =                               [K2]*[W]
%         +                              chi*[CHI]
%         +               gamma*a^3/D_eff*[K4]*[W]
%         -      A*(1-nu^2)*h^2/De*a^4/b^4*[TEMP2]
% ------------------------------------------------------------------------ 
%   D_FW_W=              [K2]+gamma*a^3/D_eff*[K4]
%         -   A*(1-nu^2)*h^2/De*a^4/b^4*D_TEMP2_DW
% ========================================================================
for p=1:MN,
TEMP1(p)=W'*C_cell{p}*W;
    for l=1:1:MN,
        D_TEMP1_DW(p,l)=C_cell{p}(l,:)*W+W'*C_cell{p}(:,l);
    end
end
TEMP1=K1\TEMP1; 
D_TEMP1_DW=K1\D_TEMP1_DW;
% ------------------------------------------------------------------------
for p=1:1:MN,
TEMP2=TEMP1'*S_cell{p}*W;
    for l=1:1:MN,
        D_TEMP2_DW(p,l)=D_TEMP1_DW(:,l)'*S_cell{p}*W+TEMP1'*S_cell{p}(:,l);
    end
end 
%  c4    = gamma*a^3/D_eff;
%  c_phi = A*(1-nu^2)*h^2/D_eff*a^4/b^4;
FW     = ( K2 + c4*K4 )*W  + chi*CHI - c_phi*TEMP2;
D_FW_W =   K2 + c4*K4      - c_phi*D_TEMP2_DW;
FW=sparse(FW);
D_FW_W=sparse(D_FW_W);