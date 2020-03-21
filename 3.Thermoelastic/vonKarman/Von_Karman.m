
% ========================================================================
% 
%   This function caculates the stiffness matrix of weak form of the von 
%   Karman plate equation. 
%  
%   The boudary conditions are simple supported for the out-of-plan
%   displacement and homogeneous for the airy stress function.
%   
%           2016 - 10 -22 by Mingguan YANG
% =======================================================================

% =======================================================================
% 
%  The weak formulation of the equation is as follows:
% 
%               This is the adimenstional form(1)  
% 
%                   W=w/h;    PHI=phi/NT/a/b;  
% 
%   int((lambda^4*PHI_4X+2*lambda^2*PHI_2X2y+PHI_4y)*PSI*dX*dY)
%   =int(A*(1-nu^2)/NT/b^3*(W_XY^2-W_XX^2*W_YY^2)*PSI*dX*dY))
%                               (2)
%   int((lambda^4*W_4X+2*lambda^2*W_2X2Y+W_4Y)*V*dX*dY)
%   -gamma*a^3/D_eff*int(((X-1)*W_X)_X*V*dX*dY)
%   =NT*a^3/D_eff/b*int((PHI_XX*W_YY+PHI_YY*W_XX-2*PHI_XY*W_XY)*V*dX*dY)
% 
% ========================================================================
%               
%     The associated Matrix form of the weak formulation above is:
% ------------------------------------------------------------------------
%                           First Equation
% ------------------------------------------------------------------------
%                                              {[W]T*[C_11]*[W]}
%                                              {       .       }
%                                              {       .       }
%           [K1]*[PHI]=A*(1-nu^2)*a*h^2/NT/b^3*{[W]T*[C_ij]*[W]}
%                                              {       .       }
%                                              {       .       }
%                                              {[W]T*[C_MN]*[W]}
% ------------------------------------------------------------------------
%                           Second Equation
% ------------------------------------------------------------------------
%                                                            {[PHI]T*[S_11]*[W]}
%                                                            {         .       }
%                                                            {         .       }
% [K2]*[W]+chi*[CHI]+gamma*a^3/D_eff*[K4]*[W]=NT*a^3/D_eff/b*{[PHI]T*[S_ij]*[W]}
%                                                            {         .       }
%                                                            {         .       }
%                                                            {[PHI]T*[S_MN]*[W]}

% ========================================================================
% ========================================================================
%                            Input parameters                   
% ========================================================================

function [K1,K2,K4,C_cell,S_cell,CHI]=Von_Karman(M,N,a,b,h)
disp('Sub Step Start: Calculation of stiffness matrix')

MN=M*N;
lambda=a/b;
K1=zeros(MN);
K2=zeros(MN);
K4=zeros(MN);

C_cell=cell(MN,1);
S_cell=cell(MN,1);
CHI=zeros(MN,1);


Ix_cc2c=cell(MN,1);
Ix_ss2c=cell(MN,1);


Iy_cc2c=cell(MN,1);
Iy_ss2c=cell(MN,1);

Ix_2css=cell(MN,1);
Ix_2scs=cell(MN,1);

Iy_2css=cell(MN,1);
Iy_2scs=cell(MN,1);
% ========================================================================
%              Creating indexing matrixes
% ========================================================================
M_row_matrix=repmat((1:1:M)',1,M);
M_column_matrix=repmat(1:1:M,M,1);
N_row_matrix=repmat((1:1:N)',1,N);
N_column_matrix=repmat(1:1:N,N,1);
disp('Step 1. Preparation')
% ========================================================================
%              1.Calculation of integral matrix in x direction
% 
%           Ix_cc2c=int(cos(m*pi*x)*cos(s*pi*x)*cos(2*i*pi*x)*dx);
%           Ix_ss2c=int(sin(m*pi*x)*sin(s*pi*x)*cos(2*i*pi*x)*dx);
% 
% ========================================================================

for i=1:1:M/2-1, 
%     Caculation of line AB(m=s+2*i) in the m-s matrix
    Ix_AB = sparse(...
                    (2*i+1):M,...
                    1:(M-2*i),...
                    1/4*ones(M-2*i,1),M,M...
                   );
%     Caculation of line CD(m+s=2*i) in the m-s matrix   
    Ix_CD = sparse(...
                    1:(2*i-1),...
                    (2*i-1):-1:1,...
                    1/4*ones(2*i-1,1),M,M...
                   );
    Ix_ss2c_matrix   = Ix_AB+Ix_AB'-Ix_CD;
    Ix_cc2c_matrix   = Ix_AB+Ix_AB'+Ix_CD;
    Ix_ss2c{i} = Ix_ss2c_matrix;
    Ix_cc2c{i} = Ix_cc2c_matrix;
end
%     Caculation of line CD(m+s=2*i) for i=M/2 in the m-s matrix  
Ix_ss2c{M/2} = sparse(...
                        1:(M-1),...
                        (M-1):-1:1,...
                        -1/4*ones(M-1,1),M,M...
                     );
Ix_cc2c{M/2}= -Ix_ss2c{M/2};
for i=M/2+1:1:M,
%     Caculation of line EF(m=s+2*i) in the m-s matrix
    Ix_EF = sparse(...
                    (2*i-M):M,...
                    M:-1:(2*i-M),...
                    1/4*ones(2*M-2*i+1,1),M,M...
                   );
    Ix_ss2c_matrix = -Ix_EF;
    Ix_cc2c_matrix = Ix_EF;
    Ix_ss2c{i} = Ix_ss2c_matrix;
    Ix_cc2c{i} = Ix_cc2c_matrix;
end

% ========================================================================
%              2.Calculation of integral matrix in y direction
% 
%           Iy_cc2c=int(cos(n*pi*y)*cos(t*pi*y)*cos(2*j*pi*y)*dy);
%           Iy_ss2c=int(sin(2*n*pi*y)*dy)*sin(n*pi*y)*cos(2*j*pi*y);
% 
% ========================================================================

for i=1:1:N/2-1, 
%     Caculation of line AB(n=t+2*j) in the n-t matrix
    Iy_AB = sparse((2*i+1):N,1:(N-2*i),1/4*ones(N-2*i,1),N,N);
%     Caculation of line CD(n+t=2*j) in the n-t matrix   
    Iy_CD = sparse(1:(2*i-1),(2*i-1):-1:1,1/4*ones(2*i-1,1),N,N);
    
    Iy_ss2c_matrix   = Iy_AB+Iy_AB'-Iy_CD;
    Iy_cc2c_matrix   = Iy_AB+Iy_AB'+Iy_CD;
    
    Iy_ss2c{i} = Iy_ss2c_matrix;
    Iy_cc2c{i} = Iy_cc2c_matrix;
end
%     Caculation of line CD(n+t=2*j) for i=M/2 in the m-s matrix  
Iy_ss2c{N/2}= sparse(1:(N-1),(N-1):-1:1,-1/4*ones(N-1,1),N,N);
Iy_cc2c{N/2}= -Iy_ss2c{N/2};
for i=N/2+1:1:N,
%     Caculation of line EF(n=t+2*j) in the n-t matrix
    Iy_EF = sparse((2*i-N):N,N:-1:(2*i-N),1/4*ones(2*N-2*i+1,1),N,N);
    Iy_ss2c_matrix = -Iy_EF;
    Iy_cc2c_matrix = Iy_EF;
    Iy_ss2c{i} = Iy_ss2c_matrix;
    Iy_cc2c{i} = Iy_cc2c_matrix;
end
% ========================================================================
%              3.Calculation of integral matrix in x direction
% 
%           Ix_2css=int(cos(2*m*pi*x)*sin(s*pi*x)*cos(i*pi*x)*dx);
%           Ix_2scs=int(sin(2*m*pi*x)*cos(s*pi*x)*cos(i*pi*x)*dx);
% 
% ========================================================================
for i=3:1:M-2, 
    flri=floor(i/2);
%     Caculation of line AB(2*m=s+i) in the m-s matrix
    Ix_AB = sparse(...
                    flri+1 : flri+M/2,...
                    2*flri-i+2 :2: 2*flri-i+M,...
                    1/4*ones(M/2,1),M,M...
                    );  
%     Caculation of line CD(2*m=i-s) in the m-s matrix   
    Ix_CD = sparse(...
                    i-flri-1:-1:1,...
                    2*flri-i+2:2:i-2,...
                    1/4*ones(i-flri-1,1),M,M...
                    );
%     Caculation of line EF(2*m=s-i) in the m-s matrix 
    Ix_EF = sparse(...
                    1:M/2+flri-i,...
                    i+2:2:M+2*flri-i,...
                    1/4*ones(M/2+flri-i,1),M,M...
                    );    
    Ix_2css{i} = -Ix_AB+Ix_CD+Ix_EF;
    Ix_2scs{i} = Ix_AB+Ix_CD-Ix_EF;
end
for i=1:2,
    flri=floor(i/2);
%     Caculation of line AB(2*m=s+i) in the m-s matrix
    Ix_AB = sparse(...
                    flri+1 : flri+M/2,...
                    2*flri-i+2 :2: 2*flri-i+M,...
                    1/4*ones(M/2,1),M,M...
                    );
%     Caculation of line EF(2*m=s-i) in the m-s matrix 
    Ix_EF = sparse(...
                    1:M/2+flri-i,...
                    i+2:2:M+2*flri-i,...
                    1/4*ones(M/2+flri-i,1),M,M...
                    ); 
    Ix_2css{i} = -Ix_AB+Ix_EF;
    Ix_2scs{i} = Ix_AB-Ix_EF;
end
for i=M-1:M,
        flri=floor(i/2);
%     Caculation of line AB(2*m=s+i) in the m-s matrix
    Ix_AB = sparse(...
                    flri+1 : flri+M/2,...
                    2*flri-i+2 :2: 2*flri-i+M,...
                    1/4*ones(M/2,1),M,M...
                    );                
%     Caculation of line CD(2*m=i-s) in the m-s matrix   
    Ix_CD = sparse(...
                    i-flri-1:-1:1,...
                    2*flri-i+2:2:i-2,...
                    1/4*ones(i-flri-1,1),M,M...
                    );
    Ix_2css{i} = -Ix_AB+Ix_CD;
    Ix_2scs{i} = Ix_AB+Ix_CD;
end
% ========================================================================
%              4.Calculation of integral matrix in y direction
% 
%           Iy_2css=int(cos(2*m*pi*y)*sin(s*pi*y)*cos(i*pi*y)*dy);
%           Iy_2scs=int(sin(2*m*pi*y)*cos(s*pi*y)*cos(i*pi*y)*dy);
% 
% ========================================================================
for i=3:1:N-2, 
    flri=floor(i/2);
%     Caculation of line AB(2*n=t+j) in the n-t matrix
    Iy_AB = sparse(...
                    flri+1 : flri+N/2,...
                    2*flri-i+2 :2: 2*flri-i+N,...
                    1/4*ones(N/2,1),N,N...
                    );                
%     Caculation of line CD(2*n=j-t) in the m-t matrix   
    Iy_CD = sparse(...
                    i-flri-1:-1:1,...
                    2*flri-i+2:2:i-2,...
                    1/4*ones(i-flri-1,1),N,N...
                    );
%     Caculation of line EF(2*n=t-j) in the n-t matrix 
    Iy_EF = sparse(...
                    1:N/2+flri-i,...
                    i+2:2:N+2*flri-i,...
                    1/4*ones(N/2+flri-i,1),N,N...
                    );    
    Iy_2css{i} = -Iy_AB+Iy_CD+Iy_EF;
    Iy_2scs{i} = Iy_AB+Iy_CD-Iy_EF;
end
for i=1:2,
    flri=floor(i/2);
%     Caculation of line AB(2*n=t+j) in the n-t matrix
    Iy_AB = sparse(...
                    flri+1 : flri+N/2,...
                    2*flri-i+2 :2: 2*flri-i+N,...
                    1/4*ones(N/2,1),N,N...
                    );
%     Caculation of line EF(2*n=t-j) in the n-t matrix 
    Iy_EF = sparse(...
                    1:N/2+flri-i,...
                    i+2:2:N+2*flri-i,...
                    1/4*ones(N/2+flri-i,1),N,N...
                    ); 
    Iy_2css{i} = -Iy_AB+Iy_EF;
    Iy_2scs{i} = Iy_AB-Iy_EF;
end
for i=N-1:N,
        flri=floor(i/2);
%     Caculation of line AB(2*n=t+j) in the n-t matrix
    Iy_AB = sparse(...
                    flri+1 : flri+N/2,...
                    2*flri-i+2 :2: 2*flri-i+N,...
                    1/4*ones(N/2,1),N,N...
                    );                
%     Caculation of line CD(2*n=j-t) in the n-t matrix   
    Iy_CD = sparse(...
                    i-flri-1:-1:1,...
                    2*flri-i+2:2:i-2,...
                    1/4*ones(i-flri-1,1),N,N...
                    );
    Iy_2css{i} = -Iy_AB+Iy_CD;
    Iy_2scs{i} = Iy_AB+Iy_CD;
end
% ========================================================================
%              5.Calculation of stiffness matrix for the equation
% 
%           Iy_2css=int(cos(2*m*pi*y)*sin(s*pi*y)*cos(i*pi*y)*dy);
%           Iy_2scs=int(sin(2*m*pi*y)*cos(s*pi*y)*cos(i*pi*y)*dy);
% 
% ========================================================================
disp('Step 2. Calculation')
for p=1:1:MN,
%     The mapping from ij to p
    j=ceil(p/M);
    i=p-(j-1)*M;
% ========================================================================
%               The left-hand side of first equation
%               int(nabla_phi*nabla_psi*dx*dy)=[K1]*[phi]
% ========================================================================
    for t=1:1:N,
            l=(t-1)*M+i;
            K1(p,l)=8*(i*pi)^4;
    end
    for s=1:1:M,
        if s~=i,
            l=(j-1)*M+s;
            K1(p,l)=8*(j*pi)^4*lambda^4;
            K4(p,l)=s*i/2*((-1)^(s+i)-1)*(s^2+i^2)/(s^2-i^2)^2;
        end
    end
    K1(p,p) = 12*(i*pi)^4+12*(j*pi*lambda)^4+8*i^2*j^2*lambda^2*pi^4;
    K2(p,p) = pi^4/4*(i^2+lambda^2*j^2)^2;
    CHI(p)=a^2/h*((-1)^i-1)*((-1)^j-1)*(lambda^2*j/i+i/j);
    K4(p,p)=-i^2*pi^2/8;
% % ========================================================================
% %               The right-hand side of first equation
% %           int(W_XY^2-W_XX^2*W_YY^2)*PSI*dX*dY))
% %               =[[w]T*[C11]*[w],...,[w]T*[Cij]*[w],...,[w]T*[CMN]*[w]]T
% % ========================================================================
%   alpha=(n-1)*M+m;
%   beta=(t-1)*M+s;
    I1 = 1/2*speye(M)-Ix_cc2c{i};
    I2 = 1/2*speye(N)-Iy_cc2c{j};
    I3 = 1/2*speye(M)-Ix_ss2c{i};
    I4 = 1/2*speye(N)-Iy_ss2c{j};
    %   C_ij(alpha,beta)=pi^4*(m*s*I1*n*t*I2-m^2*I3*t^2*I4);  
    C_cell{p} = pi^4*(...
                          kron( N_row_matrix.*N_column_matrix.*I2,...
                                M_row_matrix.*M_column_matrix.*I1 )...
                          -kron(N_column_matrix.^2.*I4,...
                                M_row_matrix.^2.*I3 )...
                     );
    I7  = Ix_2css{i};
    I10 = Iy_2css{j};
    I8  = sparse(1:N,j,1/2*ones(N,1),N,N)-Iy_2css{j};
    I9  = sparse(1:M,i,1/2*ones(M,1),M,M)-Ix_2css{i};
    I11 = Ix_2scs{i};
    I12 = Iy_2scs{j};
%     S_ij(alpha,beta)=-4*pi^4(m^2*I7*t^2*I8+s^2*I9*n^2*I10+2*m*s*I11*n*t*I12);
%   m       M_row_matrix
%   s       M_column_matrix
%   n       N_row_matrix
%   t       N_column_matrix
    S_cell{p} =-4*pi^4*(...
                           kron(N_column_matrix.^2.*I8,M_row_matrix.^2.*I7)...
                       +   kron(N_row_matrix.^2.*I10,M_column_matrix.^2.*I9)...
                       + 2*kron(...
                                N_row_matrix.*N_column_matrix.*I12,...
                                M_row_matrix.*M_column_matrix.*I11...
                                )...
                        );
end
disp('Step 3. Storage')
K1=sparse(K1);
K2=sparse(K2);
K4=sparse(K4);
CHI=sparse(CHI);
 
    