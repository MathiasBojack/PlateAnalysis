function [prob] = building_mosek(crit,thick,F,BC,B,Disc,edg_lgth,normal_loc,mesh,Uimp)
%BUILDING_MOSEK Construction of the Mosek optimization problems
if (nargin < 10)
    Uimp.A = [];
    Uimp.b = [];
end

mesh.thick_list  = zeros(mesh.NNE,1);
for j=1:length(thick)
    mesh.thick_list(mesh.el_tag==j) = thick(j);
end

nsz = 9;


NNE = mesh.NNE;
NED = mesh.NED;
nged = mesh.nged;
area = el_area(1:NNE,mesh);
sd = 15;
switch mesh.el_type
    case 'continuous'
        disc = 0;
        strain = 1;
    case 'discontinuous-only'
        disc = 1;
        strain = 0;
    case 'discontinuous'
        disc = 1;
        strain = 1;
end


[ad,wd] = quadrature_1D(nged,'trapeze');
Nedg_lin = [(1-ad')./2 (1+ad')./2];
Nedg_quad = [ad'.*(ad'-1)./2 ad'.*(1+ad')./2  1-ad'.^2];

switch crit.type
    case 'shear-stab'
        q1 = 13; % number of conic variables for Prm_strain/point
        q2 = 7; % number of conic variables for Prm_disc/point
        p1 = 3*NNE; % total number of integration points for Prm_strain        
        Nd = mesh.Nu;
        
        omega = 1;
        Tpen = 1;
        
        Cb = sparse([2 1 0;0 sqrt(3) 0;0 0 1]);
        Cs = speye(2);
        Cn = sparse(diag([2,1]));
        
        p2 = 3*NED; % total number of intergration points for Prm_disc

        c_disc = repmat(edg_lgth'./4,3,1);
        Crit_disc = [kron(sparse(diag([1,1,2])),[blkdiag(1/sqrt(3).*Cn,4/sqrt(3));sparse(2,3)]) ...
                        kron(sparse([1 0;0 1;1 1]),[sparse(3,2);thick./sqrt(3).*Cn]) sparse(15,2)];
        
%         c_disc2 = kron([5/9;8/9;5/9],edg_lgth'./2);
%         Crit_disc = [kron(sparse(diag([1,1,2])),[blkdiag(1/sqrt(3)/thick.*Cn,4/sqrt(3)/thick);sparse(2,3)]) ...
%                         kron(sparse([(1+1/sqrt(3))/2 (1-1/sqrt(3))/2;(1-1/sqrt(3))/2 (1+1/sqrt(3))/2;1/2 1/2]),[sparse(3,2);1./sqrt(3).*Cn]) sparse(15,2)];
         
        Crit_strain = [repmat([thick./sqrt(3).*Cb;sparse(8,3)],3,1) kron(speye(3),[sparse(3,3);1./sqrt(3).*Cb;sparse(5,3)]) kron(speye(3),[sparse(6,2);Cs;sparse(3,2)]) kron(speye(3),[sparse(8,1);Tpen])];

        
        ll = setdiff((1:sd*NED)',BC);
        C_disc = kron(speye(NED),Crit_disc);

                
                Iaux = [sparse(5,1) [speye(2) sparse(2,4); sparse(1,5) omega;sparse(2,2) diag([1 1]) sparse(2,2)]];
                Iaux2 = [sparse(5,1) [speye(2) sparse(2,4); sparse(1,5) omega;sparse(2,2) diag([0 0]) sparse(2,2)]];
                I1 = kron(disc*speye(NED),sparse(blkdiag(Iaux,Iaux,Iaux)));
                
                
               % I1 = kron(speye(p2),[sparse(5,1) [speye(2) sparse(2,4); sparse(1,5) 1;sparse(2,2) diag([0 0]) sparse(2,2)]]);
                
                
                M_disc = [C_disc(:,ll)*Disc(ll,:) sparse(p2*(q2-2),p1*q1) -I1];
                
                f = sparse(Nd+q1*p1+q2*p2,1);
                f(Nd+(1:q1:q1*p1)) = repmat((crit.strength.*area)'./3,3,1);
                f(Nd+(8:q1:q1*p1)) = 1e3*repmat((crit.strength.*area)'./3,3,1);
                f(Nd+q1*p1+(1:q2:q2*p2)) = c_disc;
                f(Nd+q1*p1+(6:q2:q2*p2)) = c_disc;

                % rhs members
                beq = sparse(1+(q1-2)*p1+(q2-2)*p2+5*NED,1);
                beq(1) = 1;

%                 prob.a = [F' sparse(1,q1*p1+q2*p2);...
%                                   kron(speye(NNE),Crit_strain)*B -kron(speye(p1),sparse(blkdiag([sparse(6,1) speye(6)],[sparse(2,1) speye(2)]))) sparse(p1*(q1-2),p2*q2); ...
%                                   M_disc];                
                              
               prob.a = [F' sparse(1,q1*p1+q2*p2);...
                                  kron(speye(NNE),Crit_strain)*B -kron(speye(p1),sparse(blkdiag([sparse(6,1) speye(6)],[sparse(5,1) omega*speye(5)]))) sparse(p1*(q1-2),p2*q2); ...
                                  M_disc;...
                                  [Disc(3:15:end,:);Disc(6:15:end,:);Disc(9:15:end,:);0*Tpen.*[Disc(14:15:end,:);Disc(15:15:end,:)]] sparse(5*NED,q1*p1+q2*p2)];

                prob.cones = cell(2*p1 + 2*p2,1);
                for t=1:p1
                    prob.cones{t}.type = 'MSK_CT_QUAD';
                    prob.cones{t}.sub = Nd+q1*(t-1)+(1:7);
                    prob.cones{p1+t}.type = 'MSK_CT_QUAD';
                    prob.cones{p1+t}.sub = Nd+q1*(t-1)+(8:q1);
                end

                for t2=1:p2
                    prob.cones{2*p1+t2}.type = 'MSK_CT_QUAD';
                    prob.cones{2*p1+t2}.sub = Nd+p1*q1+q2*(t2-1)+(1:5);
                    prob.cones{2*p1+p2+t2}.type = 'MSK_CT_QUAD';
                    prob.cones{2*p1+p2+t2}.sub = Nd+p1*q1+q2*(t2-1)+(6:q2);
                end
                      
                
    case 'vonMises'
        if not(isfield(crit,'red_fact'))
            kred = ones(1,crit.ngz);
        else
            kred = crit.red_fact;
        end
        
        if not(isfield(crit,'strength_ed'))
            crit.strength_ed = ones(mesh.NED,1);
        end

        
        ngz = crit.ngz;
        [az,wz] = quadrature_1D(ngz,'trapeze');

        [ad,wd] = quadrature_1D(nged,'trapeze');
        Nedg_lin = [(1-ad')./2 (1+ad')./2];
        Nedg_quad = [ad'.*(ad'-1)./2 ad'.*(1+ad')./2  1-ad'.^2];
        
        q1 = 4; % number of conic variables for Prm_strain/point
        q2 = 3; % number of conic variables for Prm_disc/point
        p1 = ngz*3*NNE; % total number of integration points for Prm_strain        
        Nd = mesh.Nu;
        
        
        Cb = sparse([2 1 0;0 sqrt(3) 0;0 0 1]);
        Cn = sparse(diag([2,1]));
        
%         Cb = sparse([2.189 0.7981 0;0 1.6903 0;0 0 1]);
%         Cn = [1 0 0;0 0 1]*Cb*[1 0;0 0;0 1];
        
        p2 = ngz*nged*NED; % total number of intergration points for Prm_disc

        c_disc = kron(wd',kron(kred'.*wz',crit.strength_ed'.*edg_lgth'./2));
        
        
        Cn_r = kron(-thick^2/4.*az',Cn./sqrt(3));
        Cn_u = kron(ones(ngz,1),[thick/2*Cn./sqrt(3),sparse(2,1)]);
        Crit_disc = [kron(Nedg_quad,Cn_u) ...
                        kron(Nedg_lin,Cn_r) sparse(ngz*nged*2,2)];
        
        Crit_c = kron(-thick^2/4.*az',Cb./sqrt(3));
        Crit_e = kron(ones(ngz,1),thick/2*Cb./sqrt(3));
                    
        Crit_strain = [repmat(Crit_c,3,1) kron(speye(3),Crit_e) sparse(9*ngz,nsz)];
        
        LK = [sparse(nsz,12) speye(nsz)];

        
            ll = setdiff((1:sd*NED)',BC);
            C_disc = kron(speye(NED),Crit_disc);
            
            Iaux = [sparse(q2-1,1) speye(q2-1)];
            I1 = kron(disc*speye(p2),Iaux);
            Iw = kron(speye(NED),[kron(speye(3),[sparse(1,2) 1]) sparse(3,6)]);
            
            M_disc = [C_disc(:,ll)*Disc(ll,:) sparse(p2*(q2-1),p1*q1) -I1];

            f = sparse(Nd+q1*p1+q2*p2,1);
            f(Nd+(1:q1:q1*p1)) = repmat(kron(kred'.*wz',(crit.strength.*area)'./3),3,1);
            f(Nd+q1*p1+(1:q2:q2*p2)) = c_disc;

            % rhs members
            beq = sparse(1+(q1-1)*p1+nsz*NNE+(q2-1)*p2+5*NED,1);
            beq(1) = 1;

           prob.a = [F' sparse(1,q1*p1+q2*p2);...
                              kron(speye(NNE),Crit_strain)*B -kron(strain*speye(p1),[sparse(q1-1,1) speye(q1-1)]) sparse(p1*(q1-1),p2*q2); ...
                              kron(speye(NNE),LK)*B sparse(nsz*NNE,q1*p1+q2*p2);...
                              M_disc;...
                              [Iw(:,ll)*Disc(ll,:);Disc(14:15:end,:);Disc(15:15:end,:)] sparse(5*NED,q1*p1+q2*p2)];

            prob.cones = cell(p1 + p2,1);
            for t=1:p1
                prob.cones{t}.type = 'MSK_CT_QUAD';
                prob.cones{t}.sub = Nd+q1*(t-1)+(1:q1);
            end

            for t2=1:p2
                prob.cones{p1+t2}.type = 'MSK_CT_QUAD';
                prob.cones{p1+t2}.sub = Nd+p1*q1+q2*(t2-1)+(1:q2);
            end
        
            
                
    case 'vonMises-varthick'
        if not(isfield(crit,'red_fact'))
            kred = ones(1,crit.ngz);
        else
            kred = crit.red_fact;
        end
        
        if not(isfield(crit,'strength_ed'))
            crit.strength_ed = ones(mesh.NED,1);
        end

        
        ngz = crit.ngz;
        [az,wz] = quadrature_1D(ngz,'trapeze');

        [ad,wd] = quadrature_1D(nged,'trapeze');
        Nedg_lin = [(1-ad')./2 (1+ad')./2];
        Nedg_quad = [ad'.*(ad'-1)./2 ad'.*(1+ad')./2  1-ad'.^2];
        
        q1 = 4; % number of conic variables for Prm_strain/point
        q2 = 3; % number of conic variables for Prm_disc/point
        p1 = ngz*3*NNE; % total number of integration points for Prm_strain        
        Nd = mesh.Nu;
        
        
        Cb = sparse([2 1 0;0 sqrt(3) 0;0 0 1]);
        Cn = sparse(diag([2,1]));
        
%         Cb = sparse([2.189 0.7981 0;0 1.6903 0;0 0 1]);
%         Cn = [1 0 0;0 0 1]*Cb*[1 0;0 0;0 1];
        
        p2 = ngz*nged*NED; % total number of intergration points for Prm_disc

        c_disc = kron(wd',kron(kred'.*wz',crit.strength_ed'.*edg_lgth'./2));
        
        
        Cn_r = kron(-az',Cn./sqrt(3));
        Cn_u = kron(ones(ngz,1),[2/thick*Cn./sqrt(3),sparse(2,1)]);
        Crit_disc = [kron(Nedg_quad,Cn_u) ...
                        kron(Nedg_lin,Cn_r) sparse(ngz*nged*2,2)];
        
        Crit_c = kron(-az',Cb./sqrt(3));
        Crit_e = kron(ones(ngz,1),2*Cb./sqrt(3));
                    
        Crit_strain = [repmat(Crit_c,3,1) kron(speye(3),Crit_e) sparse(9*ngz,6)];
        
        Ithick_s = kron(spdiags(1./mesh.thick_list,0,NNE,NNE),[sparse(9*ngz,3) kron(speye(3),kron(ones(ngz,1),ones(3,3))) sparse(9*ngz,6)]);
        Ithick_s = Ithick_s + kron(speye(NNE),[repmat(kron(ones(ngz,1),speye(3)),3,1) sparse(9*ngz,15)]);
        
        LK = [sparse(6,12) speye(6)];

        
            ll = setdiff((1:sd*NED)',BC);
            C_disc = kron(speye(NED),Crit_disc);
            
            Iaux = [sparse(q2-1,1) speye(q2-1)];
            I1 = kron(disc*speye(p2),Iaux);
            Iw = kron(speye(NED),[kron(speye(3),[sparse(1,2) 1]) sparse(3,6)]);
            
            M_disc = [C_disc(:,ll)*Disc(ll,:) sparse(p2*(q2-1),p1*q1) -I1];

            f = sparse(Nd+q1*p1+q2*p2,1);
            f(Nd+(1:q1:q1*p1)) = repmat(kron(kred'.*wz',(crit.strength.*area)'./3),3,1);
            f(Nd+q1*p1+(1:q2:q2*p2)) = c_disc;

            % rhs members
            beq = sparse(1+(q1-1)*p1+6*NNE+(q2-1)*p2+5*NED,1);
            beq(1) = 1;

           prob.a = [F' sparse(1,q1*p1+q2*p2);...
                              (kron(speye(NNE),Crit_strain).*Ithick_s)*B -kron(strain*speye(p1),[sparse(q1-1,1) speye(q1-1)]) sparse(p1*(q1-1),p2*q2); ...
                              kron(speye(NNE),LK)*B sparse(6*NNE,q1*p1+q2*p2);...
                              M_disc;...
                              [Iw(:,ll)*Disc(ll,:);Disc(14:15:end,:);Disc(15:15:end,:)] sparse(5*NED,q1*p1+q2*p2)];

            prob.cones = cell(p1 + p2,1);
            for t=1:p1
                prob.cones{t}.type = 'MSK_CT_QUAD';
                prob.cones{t}.sub = Nd+q1*(t-1)+(1:q1);
            end

            for t2=1:p2
                prob.cones{p1+t2}.type = 'MSK_CT_QUAD';
                prob.cones{p1+t2}.sub = Nd+p1*q1+q2*(t2-1)+(1:q2);
            end
        
            
               
    case 'membrane'
        if not(isfield(crit,'red_fact'))
            kred = ones(1,crit.ngz);
        else
            kred = crit.red_fact;
        end
        

        
        ngz = crit.ngz;
        [az,wz] = quadrature_1D(ngz,'trapeze');

        
        q1 = 4; % number of conic variables for Prm_strain/point
        q2 = 3; % number of conic variables for Prm_disc/point
        p1 = ngz*3*NNE; % total number of integration points for Prm_strain        
        Nd = mesh.Nu;
        
        
        Cb = sparse([2 1 0;0 sqrt(3) 0;0 0 1]);
        Cn = sparse(diag([2,1]));
        
        p2 = ngz*3*NED; % total number of intergration points for Prm_disc

        c_disc = repmat(kron(kred'.*wz',crit.strength_ed'.*edg_lgth'./4),3,1);
        
        
        Cn_r = kron(-az',thick./2.*Cn./sqrt(3));
        Cn_u = kron(ones(ngz,1),[Cn./sqrt(3),sparse(2,1)]);
        Crit_disc = [kron(sparse(diag([1,1,2])),Cn_u) ...
                        kron(sparse([1 0;0 1;1 1]),Cn_r) sparse(ngz*6,2)];
        
        Crit_c = kron(-az',thick/2.*Cb./sqrt(3));
        Crit_e = kron(ones(ngz,1),Cb./sqrt(3));
                    
        Crit_strain = [repmat(Crit_c,3,1) kron(speye(3),Crit_e) sparse(9*ngz,6)];
        LK = [sparse(6,12) speye(6)];

        
            ll = setdiff((1:sd*NED)',BC);
            C_disc = kron(speye(NED),Crit_disc);
            
            Iaux = [sparse(q2-1,1) speye(q2-1)];
            I1 = kron(disc*speye(ngz*NED),sparse(blkdiag(Iaux,Iaux,Iaux)));
            Iw = kron(speye(NED),[kron(speye(3),[sparse(1,2) 1]) sparse(3,6)]);
            
            M_disc = [C_disc(:,ll)*Disc(ll,:) sparse(p2*(q2-1),p1*q1) -I1];

            f = sparse(Nd+q1*p1+q2*p2,1);
            f(Nd+(1:q1:q1*p1)) = repmat(kron(kred'.*wz',(crit.strength.*area)'./3),3,1);
            f(Nd+q1*p1+(1:q2:q2*p2)) = c_disc;

            % rhs members
            beq = sparse(1+(q1-1)*p1+6*NNE+(q2-1)*p2+5*NED,1);
            beq(1) = 1;

           prob.a = [F' sparse(1,q1*p1+q2*p2);...
                              kron(speye(NNE),Crit_strain)*B -kron(strain*speye(p1),[sparse(q1-1,1) speye(q1-1)]) sparse(p1*(q1-1),p2*q2); ...
                              kron(speye(NNE),LK)*B sparse(6*NNE,q1*p1+q2*p2);...
                              M_disc;...
                              [Iw(:,ll)*Disc(ll,:);Disc(14:15:end,:);Disc(15:15:end,:)] sparse(5*NED,q1*p1+q2*p2)];

            prob.cones = cell(p1 + p2,1);
            for t=1:p1
                prob.cones{t}.type = 'MSK_CT_QUAD';
                prob.cones{t}.sub = Nd+q1*(t-1)+(1:q1);
            end

            for t2=1:p2
                prob.cones{p1+t2}.type = 'MSK_CT_QUAD';
                prob.cones{p1+t2}.sub = Nd+p1*q1+q2*(t2-1)+(1:q2);
            end
        
        
            
    case {'MCt','MCt-aniso'}
        switch crit.type
            case 'MCt'
                crit.fc1 = crit.fc;
                crit.fc2 = crit.fc;
                crit.ft1 = crit.ft;
                crit.ft2 = crit.ft;
                crit.alp = 0;
        end
        if not(isfield(crit,'red_fact'))
            kred = ones(1,crit.ngz);
        else
            kred = crit.red_fact;
        end
        
        Ralp = [cos(crit.alp)^2 sin(crit.alp)^2 sin(2*crit.alp)/2;...
                  sin(crit.alp)^2 cos(crit.alp)^2 -sin(2*crit.alp)/2;...
                  -sin(2*crit.alp) sin(2*crit.alp) cos(2*crit.alp)];
        
        
        ngz = crit.ngz;
        [az,wz] = quadrature_1D(ngz,'trapeze');

        q1 = 6; % number of conic variables for Prm_strain/point
        q2 = 6; % number of conic variables for Prm_disc/point
        p1 = ngz*3*NNE; % total number of integration points for Prm_strain        
        Nd = mesh.Nu;
        
        
        p2 = ngz*nged*NED; % total number of intergration points for Prm_disc

        c_disc = kron(wd',kron(kred'.*wz',edg_lgth'./2));
        
        
        Cn = speye(2); 
        Cn_r = kron(-thick/2*az',Cn);
        Cn_u = kron(ones(ngz,1),[Cn,sparse(2,1)]);
        Crit_disc = [kron(Nedg_quad,Cn_u) ...
                        kron(Nedg_lin,Cn_r) sparse(ngz*nged*2,2)];
        
                    
        Cb = Ralp;           
        Crit_c = kron(-thick/2*az',Cb);
        Crit_e = kron(ones(ngz,1),Cb);
                    
        Crit_strain = [repmat(Crit_c,3,1) kron(speye(3),Crit_e) sparse(9*ngz,6)];
        LK = [sparse(6,12) speye(6)];
        
        norm = sparse(mesh.inplane_norm);
        tang = sparse(cross(norm',mesh.edg_norm')');
        Rotn = kron(diag(norm(:,1).^2),kron(speye(ngz*nged),[1 0;0 0;0 0])) + ...
                  kron(diag(norm(:,2).^2),kron(speye(ngz*nged),[0 0;1 0;0 0])) + ...
                  kron(diag(norm(:,1).*norm(:,2)),kron(speye(ngz*nged),[0 0;0 0;1 0])) + ...
                  kron(diag(norm(:,1).*tang(:,1)),kron(speye(ngz*nged),[0 1;0 0;0 0])) + ...
                  kron(diag(norm(:,2).*tang(:,2)),kron(speye(ngz*nged),[0 0;0 1;0 0])) + ...
                  kron(diag(norm(:,1).*tang(:,2)+norm(:,2).*tang(:,1)),kron(speye(ngz*nged),[0 0;0 0;0 1/2]));
          %  Rotn = kron(speye(NED*ngz*nged),[1 0;0 0;0 1]);
        
        
            ll = setdiff((1:sd*NED)',BC);
            C_disc = kron(speye(p2),Ralp)*Rotn*kron(speye(NED),Crit_disc);
            
            Iaux = [sparse(q2-1,1) speye(q2-1)];
            
            I1 = kron(speye(ngz*NED),sparse(blkdiag(Iaux,Iaux,Iaux)));
            Iw = kron(speye(NED),[kron(speye(3),[sparse(1,2) 1]) sparse(3,6)]);
            
%            M_disc = [C_disc(:,ll)*Disc(ll,:) sparse(p2*(q2-1),p1*q1) -I1];

            f = sparse(Nd+q1*p1+q2*p2,1);
            f(Nd+(1:q1:q1*p1)) = repmat(kron(kred'.*wz',(crit.ft1.*area)'./3),3,1);
            f(Nd+(2:q1:q1*p1)) = repmat(kron(kred'.*wz',(crit.ft2.*area)'./3),3,1);
            f(Nd+(4:q1:q1*p1)) = repmat(kron(kred'.*wz',(crit.fc1.*area)'./3),3,1);
            f(Nd+(5:q1:q1*p1)) = repmat(kron(kred'.*wz',(crit.fc2.*area)'./3),3,1);
            f(Nd+q1*p1+(1:q2:q2*p2)) = crit.ft1.*c_disc;
            f(Nd+q1*p1+(2:q2:q2*p2)) = crit.ft2.*c_disc;
            f(Nd+q1*p1+(4:q2:q2*p2)) = crit.fc1.*c_disc;
            f(Nd+q1*p1+(5:q2:q2*p2)) = crit.fc2.*c_disc;

            % rhs members
            beq = sparse(1+(q1-3)*p1+6*NNE+(q2-3)*p2+5*NED,1);
            beq(1) = 1;

           prob.a = [F' sparse(1,q1*p1+q2*p2);...
                              kron(speye(NNE),Crit_strain)*B -kron(strain*speye(p1),sparse([diag([1,1,sqrt(2)]),diag([-1,-1,sqrt(2)])])) sparse(p1*(q1-3),p2*q2); ...
                              kron(speye(NNE),LK)*B sparse(6*NNE,q1*p1+q2*p2);...
                              C_disc(:,ll)*Disc(ll,:) sparse(p2*(q2-3),p1*q1) -kron(disc*speye(p2),sparse([diag([1,1,sqrt(2)]),diag([-1,-1,sqrt(2)])])) ;...
                              [Iw(:,ll)*Disc(ll,:);Disc(14:15:end,:);Disc(15:15:end,:)] sparse(5*NED,q1*p1+q2*p2)];

            prob.cones = cell(2*p1 + 2*p2,1);
            for t=1:p1
                prob.cones{t}.type = 'MSK_CT_RQUAD';
                prob.cones{t}.sub = Nd+q1*(t-1)+(1:3);
                prob.cones{p1+t}.type = 'MSK_CT_RQUAD';
                prob.cones{p1+t}.sub = Nd+q1*(t-1)+(4:6);
            end

            for t2=1:p2
                prob.cones{2*p1+t2}.type = 'MSK_CT_RQUAD';
                prob.cones{2*p1+t2}.sub = Nd+p1*q1+q2*(t2-1)+(1:3);
                prob.cones{2*p1+p2+t2}.type = 'MSK_CT_RQUAD';
                prob.cones{2*p1+p2+t2}.sub = Nd+p1*q1+q2*(t2-1)+(4:6);
            end
           
                 
    case 'Tresca'
        
        
        ngz = crit.ngz;
        [az,wz] = quadrature_1D(ngz,'trapeze');
        

        q1 = 9; % number of conic variables for Prm_strain/point
        q2 = 9; % number of conic variables for Prm_disc/point
        p1 = ngz*3*NNE; % total number of integration points for Prm_strain        
        Nd = mesh.Nu;
        
        
        p2 = ngz*nged*NED; % total number of intergration points for Prm_disc

        c_disc = kron(wd',kron(wz',edg_lgth'./2));
        
        Cn = sparse([1 0;0 0;0 1]); 
        Cn_r = kron(-thick^2/4*az',Cn);
        Cn_u = kron(ones(ngz,1),[thick/2*Cn,sparse(3,1)]);
        Crit_disc = [kron(sparse(Nedg_quad),Cn_u) ...
                        kron(sparse(Nedg_lin),Cn_r) sparse(ngz*3*nged,2)];
        
                    
        Cb = speye(3);           
        Crit_c = kron(-thick^2/4*az',Cb);
        Crit_e = kron(ones(ngz,1),thick/2*Cb);
                    
        Crit_strain = [repmat(Crit_c,3,1) kron(speye(3),Crit_e) sparse(9*ngz,6)];
        LK = [sparse(6,12) speye(6)];

        
            ll = setdiff((1:sd*NED)',BC);
            C_disc = kron(speye(NED),Crit_disc);
            
            Iaux = [sparse(q2-1,1) speye(q2-1)];
            
            I1 = kron(speye(ngz*NED),sparse(blkdiag(Iaux,Iaux,Iaux)));
            Iw = kron(speye(NED),[kron(speye(3),[sparse(1,2) 1]) sparse(3,6)]);
            
%            M_disc = [C_disc(:,ll)*Disc(ll,:) sparse(p2*(q2-1),p1*q1) -I1];

            f = sparse(Nd+q1*p1+q2*p2,1);
            f(Nd+(1:q1:q1*p1)) = repmat(kron(wz',(area)'./3),3,1);
            f(Nd+(2:q1:q1*p1)) = repmat(kron(wz',(area)'./3),3,1);
            f(Nd+(4:q1:q1*p1)) = repmat(kron(wz',(area)'./3),3,1);
            f(Nd+(5:q1:q1*p1)) = repmat(kron(wz',(area)'./3),3,1);
            f(Nd+(7:q1:q1*p1)) = repmat(kron(wz',(area)'./3),3,1);
            f(Nd+q1*p1+(1:q2:q2*p2)) = c_disc;
            f(Nd+q1*p1+(2:q2:q2*p2)) = c_disc;
            f(Nd+q1*p1+(4:q2:q2*p2)) = c_disc;
            f(Nd+q1*p1+(5:q2:q2*p2)) = c_disc;
            f(Nd+q1*p1+(7:q2:q2*p2)) = c_disc;

            % rhs members
            beq = sparse(1+(q1-6)*p1+6*NNE+(q2-6)*p2+5*NED,1);
            beq(1) = 1;

           prob.a = [F' sparse(1,q1*p1+q2*p2);...
                              kron(speye(NNE),Crit_strain)*B -kron(strain*speye(p1),sparse([diag([1,1,sqrt(2)]),diag([-1,-1,sqrt(2)]),[sparse(3,1) [1 0;-1 0;0 2]]])) sparse(p1*(q1-6),p2*q2); ...
                              kron(speye(NNE),LK)*B sparse(6*NNE,q1*p1+q2*p2);...
                              C_disc(:,ll)*Disc(ll,:) sparse(p2*(q2-6),p1*q1) -kron(disc*speye(p2),sparse([diag([1,1,sqrt(2)]),diag([-1,-1,sqrt(2)]),[sparse(3,1) [1 0;-1 0;0 2]]])) ;...
                              [Iw(:,ll)*Disc(ll,:);Disc(14:15:end,:);Disc(15:15:end,:)] sparse(5*NED,q1*p1+q2*p2)];

            prob.cones = cell(3*p1 + 3*p2,1);
            for t=1:p1
                prob.cones{t}.type = 'MSK_CT_RQUAD';
                prob.cones{t}.sub = Nd+q1*(t-1)+(1:3);
                prob.cones{p1+t}.type = 'MSK_CT_RQUAD';
                prob.cones{p1+t}.sub = Nd+q1*(t-1)+(4:6);
                prob.cones{2*p1+t}.type = 'MSK_CT_QUAD';
                prob.cones{2*p1+t}.sub = Nd+q1*(t-1)+(7:9);
            end

            for t2=1:p2
                prob.cones{3*p1+t2}.type = 'MSK_CT_RQUAD';
                prob.cones{3*p1+t2}.sub = Nd+p1*q1+q2*(t2-1)+(1:3);
                prob.cones{3*p1+p2+t2}.type = 'MSK_CT_RQUAD';
                prob.cones{3*p1+p2+t2}.sub = Nd+p1*q1+q2*(t2-1)+(4:6);
                prob.cones{3*p1+2*p2+t2}.type = 'MSK_CT_RQUAD';
                prob.cones{3*p1+2*p2+t2}.sub = Nd+p1*q1+q2*(t2-1)+(7:9);
            end       
            
    case 'MCt-r'
        if not(isfield(crit,'red_fact'))
            kcred = ones(1,crit.ngz);
        else
            kcred = crit.red_fact;
        end
        if not(isfield(crit.ac,'red_fact'))
            kared = ones(1,crit.ac.nbac);
        else
            kared = crit.ac.red_fact;
        end
        
        
        ngz = crit.ngz;
        nga = crit.ac.nbac;
        
        [ad,wd] = quadrature_1D(nged,'trapeze');
        Nedg_lin = [(1-ad')./2 (1+ad')./2];
        Nedg_quad = [ad'.*(ad'-1)./2 ad'.*(1+ad')./2  1-ad'.^2];
        
        [az,wz] = quadrature_1D(ngz,'trapeze');
        xia = crit.ac.zs;
        
        q1 = 6; % number of conic variables for Prm_strain/point
        q2 = 6; % number of conic variables for Prm_disc/point
        p1 = ngz*3*NNE; % total number of integration points for Prm_strain        
        Nd = mesh.Nu;
        
        
        p2 = ngz*nged*NED; % total number of intergration points for Prm_disc

        p3 = 2*nga*3*NNE; % total number of integration points for Prm_strain,renf     
        p4 = 2*nga*nged*NED; % total number of integration points for Prm_disc,renf     
        
        ll = setdiff((1:sd*NED)',BC);
        
        c_disc = kron(wd',kron(kcred'.*wz',edg_lgth'./2));
        
        
        Cn = sparse([1 0;0 0;0 1]); 
        Cn_r = kron(-thick^2/4*az',Cn);
        Cn_u = kron(ones(ngz,1),[thick/2*Cn,sparse(3,1)]);
        Crit_disc = [kron(sparse(Nedg_quad),Cn_u) ...
                        kron(sparse(Nedg_lin),Cn_r) sparse(ngz*3*nged,2)];
        
                    
        Cb = speye(3);           
        Crit_c = kron(-thick^2/4*az',Cb);
        Crit_e = kron(ones(ngz,1),thick/2*Cb);
        Crit_strain = [repmat(Crit_c,3,1) kron(speye(3),Crit_e) sparse(9*ngz,6)];
        
        Dir = sparse(nga,3);
        Dir_n = sparse(nga,NED);
        Dir_t = sparse(nga,NED);
        for a=1:nga
            Dir(a,:) = [crit.ac.dir(a,1)^2 crit.ac.dir(a,2)^2 crit.ac.dir(a,1)*crit.ac.dir(a,2)];
            Dir_n(a,:) = crit.ac.dir(a,:)*normal_loc';
            Dir_t(a,:) = crit.ac.dir(a,:)*[0 -1;1 0]*normal_loc';
        end
        Dir_n = repmat(Dir_n,nged,1);
        Dir_t = repmat(Dir_t,nged,1);
        
        Crit_renf_e = repmat(1,nga,3).*Dir;
        Crit_renf_c = repmat(-thick/2*xia',1,3).*Dir;

        Crit_renf_strain = [repmat(Crit_renf_c,3,1) kron(speye(3),Crit_renf_e) sparse(3*nga,6)];
        
        Cn = sparse([1 0]); 
        Cn_u = kron(ones(nga,1),[1 0 0]);
        Crit_disc_n = [kron(sparse(Nedg_quad),Cn_u) ...
                        kron(sparse(Nedg_lin),kron(-xia'*thick/4,Cn)) sparse(nged*nga,2)];
        Ct = sparse([0 1]); 
        Ct_u = kron(ones(nga,1),[0 1 0]);
        Crit_disc_t = [kron(sparse(Nedg_quad),Ct_u) ...
                        kron(sparse(Nedg_lin),kron(-xia'*thick/4,Ct))  sparse(nged*nga,2)];
        C_disc_renf = sparse(diag(Dir_n(:).^2))*kron(speye(NED),Crit_disc_n) + sparse(diag(Dir_n(:).*Dir_t(:)))*kron(speye(NED),Crit_disc_t);
        C_disc_renf = C_disc_renf(:,ll)*Disc(ll,:);
        
%         reorder = reshape((1:3*nga)',3,nga)';
%         reorder = reorder(:);
%         reorder = repmat(3*nga.*(0:NED-1),3*nga,1) + repmat(reorder,1,NED);
%         C_disc_renf = C_disc_renf(reorder(:),ll)*Disc(ll,:);
        
        LK = [sparse(6,12) speye(6)];

            C_disc = kron(speye(NED),Crit_disc);
            
            Iw = kron(speye(NED),[kron(speye(3),[sparse(1,2) 1]) sparse(3,6)]);
            

            f = sparse(Nd+p3+p4+q1*p1+q2*p2,1);
            f(Nd+p3+p4+(1:q1:q1*p1)) = repmat(kron(kcred'.*wz',(crit.ft.*area)'./3),3,1);
            f(Nd+p3+p4+(2:q1:q1*p1)) = repmat(kron(kcred'.*wz',(crit.ft.*area)'./3),3,1);
            f(Nd+p3+p4+(4:q1:q1*p1)) = repmat(kron(kcred'.*wz',(crit.fc.*area)'./3),3,1);
            f(Nd+p3+p4+(5:q1:q1*p1)) = repmat(kron(kcred'.*wz',(crit.fc.*area)'./3),3,1);
            f(Nd+p3+p4+q1*p1+(1:q2:q2*p2)) = crit.ft.*c_disc;
            f(Nd+p3+p4+q1*p1+(2:q2:q2*p2)) = crit.ft.*c_disc;
            f(Nd+p3+p4+q1*p1+(4:q2:q2*p2)) = crit.fc.*c_disc;
            f(Nd+p3+p4+q1*p1+(5:q2:q2*p2)) = crit.fc.*c_disc;
            f(Nd+(1:p3)) = crit.ac.Ns.*repmat(kron(kron(kared',[1;0]),area'./3),3,1);
            f(Nd+p3+(1:p4)) = crit.ac.Ns.*kron(wd',kron(kron(kared',[1;0]),edg_lgth'./2));
            
            % rhs members
            beq = sparse(1+(q1-3)*p1+6*NNE+(q2-3)*p2+5*NED+p3/2+p4/2,1);
            beq(1) = 1;

           prob.a = [F' sparse(1,q1*p1+q2*p2+p3+p4);...
                              kron(speye(NNE),Crit_strain)*B sparse(p1*(q1-3),p3+p4) -kron(strain*speye(p1),sparse([diag([1,1,sqrt(2)]),diag([-1,-1,sqrt(2)])])) sparse(p1*(q1-3),p2*q2); ...
                              kron(speye(NNE),LK)*B sparse(6*NNE,p3+p4+q1*p1+q2*p2);...
                              C_disc(:,ll)*Disc(ll,:) sparse(p2*(q2-3),p3+p4+p1*q1) -kron(disc*speye(p2),sparse([diag([1,1,sqrt(2)]),diag([-1,-1,sqrt(2)])])) ;...
                              [Iw(:,ll)*Disc(ll,:);Disc(14:15:end,:);Disc(15:15:end,:)] sparse(5*NED,p3+p4+q1*p1+q2*p2);...
                              kron(speye(NNE),Crit_renf_strain)*B  -kron(speye(p3/2),[0 1]) sparse(p3/2,p4+q1*p1+q2*p2)
                              C_disc_renf  sparse(p4/2,p3)  -kron(speye(p4/2),[0 1]) sparse(p4/2,q1*p1+q2*p2)];

            prob.cones = cell(2*p1 + 2*p2 + p3/2+p4/2,1);
            for t=1:p1
                prob.cones{t}.type = 'MSK_CT_RQUAD';
                prob.cones{t}.sub = Nd+p3+p4+q1*(t-1)+(1:3);
                prob.cones{p1+t}.type = 'MSK_CT_RQUAD';
                prob.cones{p1+t}.sub = Nd+p3+p4+q1*(t-1)+(4:6);
            end

            for t2=1:p2
                prob.cones{2*p1+t2}.type = 'MSK_CT_RQUAD';
                prob.cones{2*p1+t2}.sub = Nd+p3+p4+p1*q1+q2*(t2-1)+(1:3);
                prob.cones{2*p1+p2+t2}.type = 'MSK_CT_RQUAD';
                prob.cones{2*p1+p2+t2}.sub = Nd+p3+p4+p1*q1+q2*(t2-1)+(4:6);
            end
            
            for t=1:p3/2
                prob.cones{2*p1+2*p2+t}.type = 'MSK_CT_QUAD';
                prob.cones{2*p1+2*p2+t}.sub = Nd+2*(t-1)+(1:2);
            end
            for t=1:p4/2
                prob.cones{2*p1+2*p2+p3/2+t}.type = 'MSK_CT_QUAD';
                prob.cones{2*p1+2*p2+p3/2+t}.sub = Nd+p3+2*(t-1)+(1:2);
            end
           
           % prob .blc = [beq(1:1+(q1-3)*p1+6*NNE+(q2-3)*p2+5*NED);-inf*ones(2*p3+2*p4,1)];
           
            
    case 'CHell-r'
        if not(isfield(crit,'red_fact'))
            kcred = ones(1,crit.ngz);
        else
            kcred = crit.red_fact;
        end
        if not(isfield(crit.ac,'red_fact'))
            kared = ones(1,crit.ac.nbac);
        else
            kared = crit.ac.red_fact;
        end
        
        
        ngz = crit.ngz;
        nga = crit.ac.nbac;
        
        [ad,wd] = quadrature_1D(nged,'trapeze');
        Nedg_lin = [(1-ad')./2 (1+ad')./2];
        Nedg_quad = [ad'.*(ad'-1)./2 ad'.*(1+ad')./2  1-ad'.^2];
        
        [az,wz] = quadrature_1D(ngz,'trapeze');
        xia = crit.ac.zs;
        
        
        param_ell = crit.param_ell;
        nell = size(param_ell,2);
        
        q1 = 4*nell+1; % number of conic variables for Prm_strain/point
        qq1 = 3*nell;
        q2 = 4*nell+1; % number of conic variables for Prm_disc/point
        qq2 = 3*nell;
        p1 = ngz*3*NNE; % total number of integration points for Prm_strain        
        Nd = mesh.Nu;
        
        p2 = ngz*nged*NED; % total number of intergration points for Prm_disc

        p3 = 2*nga*3*NNE; % total number of integration points for Prm_strain,renf     
        p4 = 2*nga*nged*NED; % total number of integration points for Prm_disc,renf     
               
        
        [Cb,centers] = mat_C_ell(param_ell);  
        
        I1 = [kron(speye(nell),[sparse(3,1) speye(3)]) sparse(3*nell,1)];
        I2 = [kron(speye(nell),sparse([1 0 0 0])) -ones(nell,1)];
               
        
        ll = setdiff((1:sd*NED)',BC);
        
        c_disc = kron(wd',kron(kcred'.*wz',edg_lgth'./2));
        
        
        Cn = Cb*sparse([1 0;0 0;0 1]); 
        Cn_r = kron(-thick^2/4*az',Cn);
        Cn_u = kron(ones(ngz,1),[thick/2*Cn,sparse(3*nell,1)]);
        Crit_disc = [kron(sparse(Nedg_quad),Cn_u) ...
                        kron(sparse(Nedg_lin),Cn_r) sparse(ngz*3*nell*nged,2)];
                    
        cent_r = kron(-thick^2/4*az',centers'*sparse([1 0;0 0;0 1]));
        cent_u = kron(ones(ngz,1),[thick/2*centers'*sparse([1 0;0 0;0 1]),sparse(nell,1)]);
        cent_disc = [kron(sparse(Nedg_quad),cent_u) ...
                        kron(sparse(Nedg_lin),cent_r) sparse(ngz*nell*nged,2)];
        
                                
        Crit_c = kron(-thick^2/4*az',Cb);
        Crit_e = kron(ones(ngz,1),thick/2*Cb);
        Crit_strain = [repmat(Crit_c,3,1) kron(speye(3),Crit_e) sparse(9*nell*ngz,6)];
        
        cent_c = kron(-thick^2/4*az',centers');
        cent_e = kron(ones(ngz,1),thick/2*centers');
        cent_strain = [repmat(cent_c,3,1) kron(speye(3),cent_e) sparse(3*nell*ngz,6)];
        
        Dir = sparse(nga,3);
        Dir_n = sparse(nga,NED);
        Dir_t = sparse(nga,NED);
        for a=1:nga
            Dir(a,:) = [crit.ac.dir(a,1)^2 crit.ac.dir(a,2)^2 crit.ac.dir(a,1)*crit.ac.dir(a,2)];
            Dir_n(a,:) = crit.ac.dir(a,:)*normal_loc';
            Dir_t(a,:) = crit.ac.dir(a,:)*[0 -1;1 0]*normal_loc';
        end
        Dir_n = repmat(Dir_n,nged,1);
        Dir_t = repmat(Dir_t,nged,1);
        
        Crit_renf_e = repmat(1,nga,3).*Dir;
        Crit_renf_c = repmat(-thick/2*xia',1,3).*Dir;

        Crit_renf_strain = [repmat(Crit_renf_c,3,1) kron(speye(3),Crit_renf_e) sparse(3*nga,6)];
        
        Cn = sparse([1 0]); 
        Cn_u = kron(ones(nga,1),[1 0 0]);
        Crit_disc_n = [kron(sparse(Nedg_quad),Cn_u) ...
                        kron(sparse(Nedg_lin),kron(-xia'*thick/4,Cn)) sparse(nged*nga,2)];
        Ct = sparse([0 1]); 
        Ct_u = kron(ones(nga,1),[0 1 0]);
        Crit_disc_t = [kron(sparse(Nedg_quad),Ct_u) ...
                        kron(sparse(Nedg_lin),kron(-xia'*thick/4,Ct))  sparse(nged*nga,2)];
        C_disc_renf = sparse(diag(Dir_n(:).^2))*kron(speye(NED),Crit_disc_n) + sparse(diag(Dir_n(:).*Dir_t(:)))*kron(speye(NED),Crit_disc_t);
        C_disc_renf = C_disc_renf(:,ll)*Disc(ll,:);
        
        
        LK = [sparse(6,12) speye(6)];

            C_disc = kron(speye(NED),Crit_disc);
            Cent_disc = kron(speye(NED),cent_disc);
            
            Iw = kron(speye(NED),[kron(speye(3),[sparse(1,2) 1]) sparse(3,6)]);
            

            f = sparse(Nd+p3+p4+q1*p1+q2*p2,1);
            f(Nd+p3+p4+(q1:q1:q1*p1)) = repmat(kron(kcred'.*wz',(area)'./3),3,1);
            f(Nd+p3+p4+q1*p1+(q2:q2:q2*p2)) = crit.fc.*c_disc;
            f(Nd+(1:p3)) = crit.ac.Ns.*repmat(kron(kron(kared',[1;0]),area'./3),3,1);
            f(Nd+p3+(1:p4)) = crit.ac.Ns.*kron(wd',kron(kron(kared',[1;0]),edg_lgth'./2));
            
            % rhs members
            buc = sparse(1+qq1*p1+6*NNE+qq2*p2+5*NED+p3/2+p4/2+p1*nell+p2*nell,1);
            buc(1) = 1;
            blc = buc;
            blc(1+qq1*p1+6*NNE+qq2*p2+5*NED+p3/2+p4/2+(1:nell*(p1+p2))) = -inf*ones(nell*(p1+p2),1);
                

           prob.a = [F' sparse(1,q1*p1+q2*p2+p3+p4);...
                              kron(speye(NNE),Crit_strain)*B sparse(p1*qq1,p3+p4) -kron(strain*speye(p1),I1) sparse(p1*qq1,p2*q2); ...
                              kron(speye(NNE),LK)*B sparse(6*NNE,p3+p4+q1*p1+q2*p2);...
                              C_disc(:,ll)*Disc(ll,:) sparse(p2*qq2,p3+p4+p1*q1) -kron(disc*speye(p2),I1) ;...
                              [Iw(:,ll)*Disc(ll,:);Disc(14:15:end,:);Disc(15:15:end,:)] sparse(5*NED,p3+p4+q1*p1+q2*p2);...
                              kron(speye(NNE),Crit_renf_strain)*B  -kron(speye(p3/2),[0 1]) sparse(p3/2,p4+q1*p1+q2*p2)
                              C_disc_renf  sparse(p4/2,p3)  -kron(speye(p4/2),[0 1]) sparse(p4/2,q1*p1+q2*p2);...
                              kron(speye(NNE),cent_strain)*B sparse(p1*nell,p3+p4) kron(strain*speye(p1),I2) sparse(p1*nell,p2*q2);...
                              Cent_disc(:,ll)*Disc(ll,:) sparse(p2*nell,p3+p4+q1*p1) kron(disc*speye(p2),I2)];

            prob.cones = cell(nell*p1 + nell*p2 + p3/2+p4/2,1);
            for t=1:p1
                for tt=1:nell
                    prob.cones{nell*(t-1)+tt}.type = 'MSK_CT_QUAD';
                    prob.cones{nell*(t-1)+tt}.sub = Nd+p3+p4+q1*(t-1)+4*(tt-1)+(1:4);
                end
            end

            for t2=1:p2
                for tt=1:nell
                    prob.cones{nell*p1+nell*(t2-1)+tt}.type = 'MSK_CT_QUAD';
                    prob.cones{nell*p1+nell*(t2-1)+tt}.sub = Nd+p3+p4+p1*q1+q2*(t2-1)+4*(tt-1)+(1:4);
                end
            end
            
            for t=1:p3/2
                prob.cones{nell*p1+nell*p2+t}.type = 'MSK_CT_QUAD';
                prob.cones{nell*p1+nell*p2+t}.sub = Nd+2*(t-1)+(1:2);
            end
            for t=1:p4/2
                prob.cones{nell*p1+nell*p2+p3/2+t}.type = 'MSK_CT_QUAD';
                prob.cones{nell*p1+nell*p2+p3/2+t}.sub = Nd+p3+2*(t-1)+(1:2);
            end
           
    case 'CHell'
        
        [ad,wd] = quadrature_1D(nged,'trapeze');
        Nedg_lin = [(1-ad')./2 (1+ad')./2];
        Nedg_quad = [ad'.*(ad'-1)./2 ad'.*(1+ad')./2  1-ad'.^2];
        
  
        param_ell = crit.param_ell;
        nell = size(param_ell,2);
        
        q1 = 7*nell+1; % number of conic variables for Prm_strain/point
        qq1 = 6*nell;
        q2 = 4*nell+1; % number of conic variables for Prm_disc/point
        qq2 = 3*nell;
        
        p1 = 3*NNE; % total number of integration points for Prm_strain        
        Nd = mesh.Nu;
        
        p2 = nged*NED; % total number of intergration points for Prm_disc

        
        [C,centers] = mat_C_ell(param_ell); 
        idx_disc = [1,3,4];
        Cdisc = mat_Cdisc_ell(C,idx_disc);
        
        Cb = C(:,4:6);
        Cm = C(:,1:3);
        centb = centers(4:6,:)';
        centm = centers(1:3,:)';
        
        Cdiscr = Cdisc(:,3);
        Cdiscu = Cdisc(:,1:2);
        centdiscr = centers(idx_disc(3),:)';
        centdiscu = centers(idx_disc([1,2]),:)';
        
        I1 = [kron(speye(nell),[sparse(6,1) speye(6)]) sparse(6*nell,1)];
        I2 = [kron(speye(nell),[1 sparse(1,6)]) -ones(nell,1)];
        Id1 = [kron(speye(nell),[sparse(3,1) speye(3)]) sparse(3*nell,1)];
        Id2 = [kron(speye(nell),[1 sparse(1,3)]) -ones(nell,1)];
        
        
        c_disc = kron(wd',edg_lgth'./2);
        
        
        Crit_disc = [kron(sparse(Nedg_quad),[Cdiscu,sparse(3*nell,1)]) ...
                        kron(sparse(Nedg_lin),[Cdiscr,sparse(3*nell,1)]) sparse(3*nell*nged,2)];
        cent_disc = [kron(sparse(Nedg_quad),[centdiscu,sparse(nell,1)]) ...
                        kron(sparse(Nedg_lin),[centdiscr,sparse(nell,1)]) sparse(nell*nged,2)];
        
                 
        Crit_strain = [repmat(Cb,3,1) kron(speye(3),Cm) sparse(18*nell,6)];
        cent_strain = [repmat(centb,3,1) kron(speye(3),centm) sparse(3*nell,6)];
        
        LK = [sparse(6,12) speye(6)];

        
            ll = setdiff((1:sd*NED)',BC);
            C_disc = kron(speye(NED),Crit_disc);
            Cent_disc = kron(speye(NED),cent_disc);
            
            Iw = kron(speye(NED),[kron(speye(3),[sparse(1,2) 1]) sparse(3,6)]);
            

            f = sparse(Nd+q1*p1+q2*p2,1);
            f(Nd+(q1:q1:q1*p1)) = repmat(area'./3,3,1);
            f(Nd+q1*p1+(q2:q2:q2*p2)) = c_disc;

            % rhs members
            buc = sparse(1+qq1*p1+6*NNE+qq2*p2+5*NED+p1*nell+p2*nell,1);
            buc(1) = 1;
            blc = buc;
            blc(1+qq1*p1+6*NNE+qq2*p2+5*NED+(1:nell*(p1+p2))) = -inf*ones(nell*(p1+p2),1);

           prob.a = [F' sparse(1,q1*p1+q2*p2);...
                              kron(speye(NNE),Crit_strain)*B -kron(strain*speye(p1),I1) sparse(p1*qq1,p2*q2); ...
                              kron(speye(NNE),LK)*B sparse(6*NNE,q1*p1+q2*p2);...
                              C_disc(:,ll)*Disc(ll,:) sparse(p2*qq2,p1*q1) -kron(disc*speye(p2),Id1) ;...
                              [Iw(:,ll)*Disc(ll,:);Disc(14:15:end,:);Disc(15:15:end,:)] sparse(5*NED,q1*p1+q2*p2);...
                              kron(speye(NNE),cent_strain)*B kron(strain*speye(p1),I2) sparse(p1*nell,p2*q2);...
                              Cent_disc(:,ll)*Disc(ll,:) sparse(p2*nell,q1*p1) kron(disc*speye(p2),Id2)];

            prob.cones = cell(nell*p1 + nell*p2,1);
            for t=1:p1
                for tt=1:nell
                    prob.cones{nell*(t-1)+tt}.type = 'MSK_CT_QUAD';
                    prob.cones{nell*(t-1)+tt}.sub = Nd+q1*(t-1)+7*(tt-1)+(1:7);
                end
            end

            for t2=1:p2
                for tt=1:nell
                    prob.cones{nell*p1+nell*(t2-1)+tt}.type = 'MSK_CT_QUAD';
                    prob.cones{nell*p1+nell*(t2-1)+tt}.sub = Nd+p1*q1+q2*(t2-1)+4*(tt-1)+(1:4);
                end
            end
        
     case 'CHell-nodisc'
        
        
  
        param_ell = crit.param_ell;
        nell = size(param_ell,2);
        
        q1 = 7*nell+1; % number of conic variables for Prm_strain/point
        qq1 = 6*nell;
        
        p1 = 3*NNE; % total number of integration points for Prm_strain        
        Nd = mesh.Nu;
        
        p2 = nged*NED; % total number of intergration points for Prm_disc

        
        [C,centers] = mat_C_ell(param_ell); 
        Cb = C(:,4:6);
        Cm = C(:,1:3);
        centb = centers(4:6,:)';
        centm = centers(1:3,:)';
        
        I1 = [kron(speye(nell),[sparse(6,1) speye(6)]) sparse(6*nell,1)];
        I2 = [kron(speye(nell),[1 sparse(1,6)]) -ones(nell,1)];
                 
        Crit_strain = [repmat(Cb,3,1) kron(speye(3),Cm) sparse(18*nell,6)];
        cent_strain = [repmat(centb,3,1) kron(speye(3),centm) sparse(3*nell,6)];
        
        LK = [sparse(6,12) speye(6)];

        
            ll = setdiff((1:sd*NED)',BC);
            
            Iw = speye(sd*NED);
            

            f = sparse(Nd+q1*p1,1);
            f(Nd+(q1:q1:q1*p1)) = repmat(area'./3,3,1);

            % rhs members
            buc = sparse(1+qq1*p1+6*NNE+sd*NED+p1*nell,1);
            buc(1) = 1;
            blc = buc;
            blc(1+qq1*p1+6*NNE+sd*NED+(1:nell*p1)) = -inf*ones(nell*p1,1);

           prob.a = [F' sparse(1,q1*p1);...
                              kron(speye(NNE),Crit_strain)*B -kron(strain*speye(p1),I1);...
                              kron(speye(NNE),LK)*B sparse(6*NNE,q1*p1);...
                              Iw(:,ll)*Disc(ll,:) sparse(sd*NED,q1*p1);...
                              kron(speye(NNE),cent_strain)*B kron(strain*speye(p1),I2)];

            prob.cones = cell(nell*p1,1);
            for t=1:p1
                for tt=1:nell
                    prob.cones{nell*(t-1)+tt}.type = 'MSK_CT_QUAD';
                    prob.cones{nell*(t-1)+tt}.sub = Nd+q1*(t-1)+7*(tt-1)+(1:7);
                end
            end

            
    case 'MCt-r2'
        if not(isfield(crit,'red_fact'))
            kcred = ones(1,crit.ngz);
        else
            kcred = crit.red_fact;
        end
        if not(isfield(crit.ac,'red_fact'))
            kared = ones(1,crit.ac.nbac);
        else
            kared = crit.ac.red_fact;
        end
        if not(isfield(crit,'red_fact_t'))
            ktred = kcred;
        else
            ktred = crit.red_fact_t;
        end
        
        
        ngz = crit.ngz;
        nga = crit.ac.nbac;
        
        [ad,wd] = quadrature_1D(nged,'trapeze');
        Nedg_lin = [(1-ad')./2 (1+ad')./2];
        Nedg_quad = [ad'.*(ad'-1)./2 ad'.*(1+ad')./2  1-ad'.^2];
        
        [az,wz] = quadrature_1D(ngz,'trapeze');
        xia = crit.ac.zs;
        
        q1 = 6; % number of conic variables for Prm_strain/point
        q2 = 6; % number of conic variables for Prm_disc/point
        p1 = ngz*3*NNE; % total number of integration points for Prm_strain        
        Nd = mesh.Nu;
        
        
        p2 = ngz*nged*NED; % total number of intergration points for Prm_disc

        p3 = nga*3*NNE; % total number of integration points for Prm_strain,renf     
        p4 = nga*nged*NED; % total number of integration points for Prm_disc,renf     
        
        ll = setdiff((1:sd*NED)',BC);
        
        c_disc_c = kron(wd',kron(kcred'.*wz',edg_lgth'./2));
        c_disc_t = kron(wd',kron(ktred'.*wz',edg_lgth'./2));
        
        
        Cn = sparse([1 0;0 1]); 
        Cn_r = kron(-thick^2/4*az',Cn);
        Cn_u = kron(ones(ngz,1),[thick/2*Cn,sparse(2,1)]);
        Crit_disc = [kron(sparse(Nedg_quad),Cn_u) ...
                        kron(sparse(Nedg_lin),Cn_r) sparse(ngz*2*nged,2)];
        
                    
        Cb = speye(3);           
        Crit_c = kron(-thick^2/4*az',Cb);
        Crit_e = kron(ones(ngz,1),thick/2*Cb);
        Crit_strain = [repmat(Crit_c,3,1) kron(speye(3),Crit_e) sparse(9*ngz,6)];
        
        Dir = sparse(nga,3);
        Dir_n = sparse(nga,NED);
        Dir_t = sparse(nga,NED);
        for a=1:nga
            Dir(a,:) = [crit.ac.dir(a,1)^2 crit.ac.dir(a,2)^2 crit.ac.dir(a,1)*crit.ac.dir(a,2)];
            Dir_n(a,:) = crit.ac.dir(a,:)*normal_loc';
            Dir_t(a,:) = crit.ac.dir(a,:)*[0 -1;1 0]*normal_loc';
        end
        Dir_n = repmat(Dir_n,nged,1);
        Dir_t = repmat(Dir_t,nged,1);
        
        Crit_renf_e = repmat(1,nga,3).*Dir;
        Crit_renf_c = repmat(-thick/2*xia',1,3).*Dir;

        Crit_renf_strain = [repmat(Crit_renf_c,3,1) kron(speye(3),Crit_renf_e) sparse(3*nga,6)];
        
        Cn = sparse([1 0]); 
        Cn_u = kron(ones(nga,1),[1 0 0]);
        Crit_disc_n = [kron(sparse(Nedg_quad),Cn_u) ...
                        kron(sparse(Nedg_lin),kron(-xia'*thick/2,Cn)) sparse(nged*nga,2)];
        Ct = sparse([0 1]); 
        Ct_u = kron(ones(nga,1),[0 1 0]);
        Crit_disc_t = [kron(sparse(Nedg_quad),Ct_u) ...
                        kron(sparse(Nedg_lin),kron(-xia'*thick/2,Ct))  sparse(nged*nga,2)];
        C_disc_renf = sparse(diag(Dir_n(:).^2))*kron(speye(NED),Crit_disc_n) + sparse(diag(Dir_n(:).*Dir_t(:)))*kron(speye(NED),Crit_disc_t);
        C_disc_renf = C_disc_renf(:,ll)*Disc(ll,:);

        
        LK = [sparse(6,12) speye(6)];

            Rotn = kron(speye(NED*ngz*nged),[1 0;0 0;0 1]);
        
        norm = sparse(mesh.inplane_norm);
        tang = sparse(cross(norm',mesh.edg_norm')');
        Rotn = kron(diag(norm(:,1).^2),kron(speye(ngz*nged),[1 0;0 0;0 0])) + ...
                  kron(diag(norm(:,2).^2),kron(speye(ngz*nged),[0 0;1 0;0 0])) + ...
                  kron(diag(norm(:,1).*norm(:,2)),kron(speye(ngz*nged),[0 0;0 0;1 0])) + ...
                  kron(diag(norm(:,1).*tang(:,1)),kron(speye(ngz*nged),[0 1;0 0;0 0])) + ...
                  kron(diag(norm(:,2).*tang(:,2)),kron(speye(ngz*nged),[0 0;0 1;0 0])) + ...
                  kron(diag(norm(:,1).*tang(:,2)+norm(:,2).*tang(:,1)),kron(speye(ngz*nged),[0 0;0 0;0 1/2]));
              
            C_disc = Rotn*kron(speye(NED),Crit_disc);
            
            Iw = kron(speye(NED),[kron(speye(3),[sparse(1,2) 1]) sparse(3,6)]);
            

            f = sparse(Nd+p3+p4+q1*p1+q2*p2,1);
            f(Nd+p3+p4+(1:q1:q1*p1)) = repmat(kron(ktred'.*wz',(crit.ft.*area)'./3),3,1);
            f(Nd+p3+p4+(2:q1:q1*p1)) = repmat(kron(ktred'.*wz',(crit.ft.*area)'./3),3,1);
            f(Nd+p3+p4+(4:q1:q1*p1)) = repmat(kron(kcred'.*wz',(crit.fc.*area)'./3),3,1);
            f(Nd+p3+p4+(5:q1:q1*p1)) = repmat(kron(kcred'.*wz',(crit.fc.*area)'./3),3,1);
            f(Nd+p3+p4+q1*p1+(1:q2:q2*p2)) = crit.ft.*c_disc_t;
            f(Nd+p3+p4+q1*p1+(2:q2:q2*p2)) = crit.ft.*c_disc_t;
            f(Nd+p3+p4+q1*p1+(4:q2:q2*p2)) = crit.fc.*c_disc_c;
            f(Nd+p3+p4+q1*p1+(5:q2:q2*p2)) = crit.fc.*c_disc_c;
            f(Nd+(1:p3)) = repmat(kron(crit.ac.Ns'.*kared',area'./3),3,1);
            f(Nd+p3+(1:p4)) = kron(wd',kron(crit.ac.Ns'.*kared',edg_lgth'./2));
            
            % rhs members
            beq = sparse(1+(q1-3)*p1+6*NNE+(q2-3)*p2+5*NED,1);
            beq(1) = 1;

           prob.a = [F' sparse(1,q1*p1+q2*p2+p3+p4);...
                              kron(speye(NNE),Crit_strain)*B sparse(p1*(q1-3),p3+p4) -kron(strain*speye(p1),sparse([diag([1,1,sqrt(2)]),diag([-1,-1,sqrt(2)])])) sparse(p1*(q1-3),p2*q2); ...
                              kron(speye(NNE),LK)*B sparse(6*NNE,p3+p4+q1*p1+q2*p2);...
                              C_disc(:,ll)*Disc(ll,:) sparse(p2*(q2-3),p3+p4+p1*q1) -kron(disc*speye(p2),sparse([diag([1,1,sqrt(2)]),diag([-1,-1,sqrt(2)])])) ;...
                              [Iw(:,ll)*Disc(ll,:);Disc(14:15:end,:);Disc(15:15:end,:)] sparse(5*NED,p3+p4+q1*p1+q2*p2);...
                              kron(speye(NNE),kron(Crit_renf_strain,[1;-1]))*B  -kron(speye(p3),[1;1]) sparse(2*p3,p4+q1*p1+q2*p2)
                              kron(C_disc_renf,[1;-1])  sparse(2*p4,p3)  -kron(speye(p4),[1;1]) sparse(2*p4,q1*p1+q2*p2)];

            prob.cones = cell(2*p1 + 2*p2,1);
            for t=1:p1
                prob.cones{t}.type = 'MSK_CT_RQUAD';
                prob.cones{t}.sub = Nd+p3+p4+q1*(t-1)+(1:3);
                prob.cones{p1+t}.type = 'MSK_CT_RQUAD';
                prob.cones{p1+t}.sub = Nd+p3+p4+q1*(t-1)+(4:6);
            end

            for t2=1:p2
                prob.cones{2*p1+t2}.type = 'MSK_CT_RQUAD';
                prob.cones{2*p1+t2}.sub = Nd+p3+p4+p1*q1+q2*(t2-1)+(1:3);
                prob.cones{2*p1+p2+t2}.type = 'MSK_CT_RQUAD';
                prob.cones{2*p1+p2+t2}.sub = Nd+p3+p4+p1*q1+q2*(t2-1)+(4:6);
            end
            
           
           blc = [beq;-inf*ones(2*p3+2*p4,1)];
           buc = [beq;zeros(2*p3+2*p4,1)];
           
           
    case 'Ilyushin'
        
        q1 = 7; % number of conic variables for Prm_strain/point
        q2 = 5; % number of conic variables for Prm_disc/point
        p1 = 2*3*NNE; % total number of integration points for Prm_strain        
        Nd = mesh.Nu;
        
        p2 = 2*nged*NED; % total number of intergration points for Prm_disc

        [ad,wd] = quadrature_1D(nged,'trapeze');
        Nedg_lin = [(1-ad')./2 (1+ad')./2];
        Nedg_quad = [ad'.*(ad'-1)./2 ad'.*(1+ad')./2  1-ad'.^2];
        
        P = sparse([1 -1/2 0;0 sqrt(3)/2 0; 0 0 sqrt(3)]);
        TP1 = [P sqrt(3)/6.*P; sparse(3,3) sqrt(33)/6.*P]'*diag([ones(1,3)./thick,ones(1,3).*4/thick^2]);
        TP2 = [P -sqrt(3)/6.*P; sparse(3,3) sqrt(33)/6.*P]'*diag([ones(1,3)./thick,ones(1,3).*4/thick^2]);
        tp1 = TP1([1 3 4 6],[1 3 4 6]);
        tp2 = TP2([1 3 4 6],[1 3 4 6]);
        
        
        
        c_disc = kron(wd',edg_lgth'./2);
        
        
        Cn_r = [sparse(2,2);speye(2)];
        Cn_u = [speye(2) sparse(2,1);sparse(2,3)];
        Crit_disc = [kron(Nedg_quad,Cn_u) ...
                        kron(Nedg_lin,Cn_r) sparse(4*nged,2)];
                            
        Crit_strain = [repmat([sparse(3,3);speye(3)],3,1) kron(speye(3),[speye(3);sparse(3,3)]) sparse(18,6)];
        
        LK = [sparse(6,12) speye(6)];

        
            ll = setdiff((1:sd*NED)',BC);
            C_disc = kron(speye(NED),Crit_disc);
            
            Iaux = [sparse(q2-1,1) tp1 sparse(q2-1,1) tp2];
            
            I1 = kron(disc*speye(p2/2),Iaux);
            Iw = kron(speye(NED),[kron(speye(3),[sparse(1,2) 1]) sparse(3,6)]);
            
            M_disc = [C_disc(:,ll)*Disc(ll,:) sparse(p2/2*(q2-1),p1*q1) -I1];

            f = sparse(Nd+q1*p1+q2*p2,1);
            
            c_strain = repmat(area'./3,3,1);
            f(Nd+(1:q1:q1*p1)) = kron(c_strain(:),ones(2,1));
            f(Nd+q1*p1+(1:q2:q2*p2)) = kron(c_disc(:),ones(2,1));

            % rhs members
            beq = sparse(1+(q1-1)*p1/2+6*NNE+(q2-1)*p2/2+5*NED,1);
            beq(1) = 1;

           prob.a = [F' sparse(1,q1*p1+q2*p2);...
                              kron(speye(NNE),Crit_strain)*B -kron(strain*speye(p1/2),[sparse(q1-1,1) TP1 sparse(q1-1,1) TP2]) sparse(p1/2*(q1-1),p2*q2); ...
                              kron(speye(NNE),LK)*B sparse(6*NNE,q1*p1+q2*p2);...
                              M_disc;...
                              [Iw(:,ll)*Disc(ll,:);Disc(14:15:end,:);Disc(15:15:end,:)] sparse(5*NED,q1*p1+q2*p2)];

            prob.cones = cell(p1 + p2,1);
            for t=1:p1
                prob.cones{t}.type = 'MSK_CT_QUAD';
                prob.cones{t}.sub = Nd+q1*(t-1)+(1:q1);
            end

            for t2=1:p2
                prob.cones{p1+t2}.type = 'MSK_CT_QUAD';
                prob.cones{p1+t2}.sub = Nd+p1*q1+q2*(t2-1)+(1:q2);
            end
        
                       
end  

        prob.c = f;
%        prob .blx = [-inf*ones(Nd,1);sparse(p3+p4,1);-inf*ones(q1*p1+q2*p2,1)];
switch crit.type
    case {'CHell','CHell-nodisc','MCt-r2'}
        prob.blc = blc;
        prob.buc = buc;
    otherwise
        prob .blc = beq;
        prob .buc = beq;     
end
        prob .blx = [];
        prob .bux = [];  

end

