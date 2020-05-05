function prob = building_mosek(crit,thick,Ncheck,NM,NNV,Fu,Fl,H,Cont,mesh)
%BUILDING_MOSEK Summary of this function goes here
%   Detailed explanation goes here

NNE = mesh.NNE;
neq = size(H,1)+size(Cont,1);
Nu = size(H,2);


switch crit.type
        
        case 'vonMises'
        if not(isfield(crit,'red_fact'))
            kred = ones(1,crit.ngz);
        else
            kred = crit.red_fact;
        end
        
        if (size(thick,1)~=mesh.NNE)
            thick = thick.*ones(mesh.NNE);
        end
               
        ngz = crit.ngz;
        [az,wz] = quadrature_1D(ngz,'uniform');
        
        q = 4;
        p = ngz*Ncheck*NNE;
        
        
        Cb = sparse([sqrt(3) 1 0;0 2 0;0 0 1]./sqrt(3));
        AN = kron(kred.*wz,[sparse(3,1) Cb]);
        AM = kron(-kred.*wz.*az,[sparse(3,1) Cb]);
        % Crit = kron(speye(Ncheck*NNE),[AN;AM]);
        
        Crit = kron(kron(spdiags(2./thick,0,NNE,NNE),speye(Ncheck)),[AN;0.*AM])+kron(speye(Ncheck*NNE),[0.*AN;AM]);
        
        
        NN = sparse(Ncheck*6,33);
        
        NN(1:6:6*Ncheck,1:3:9) = NNV(1:Ncheck,:);
        NN(2:6:6*Ncheck,2:3:9) = NNV(1:Ncheck,:);
        NN(3:6:6*Ncheck,3:3:9) = NNV(1:Ncheck,:);
        NN(4:6:6*Ncheck,9+(1:3:18)) = NM(1:Ncheck,:);
        NN(5:6:6*Ncheck,9+(2:3:18)) = NM(1:Ncheck,:);
        NN(6:6:6*Ncheck,9+(3:3:18)) = NM(1:Ncheck,:);
        
        I_NM = kron(speye(NNE),NN);        
        

        % rhs members
        beq = [sparse(neq+2*(q-1)*Ncheck*NNE,1);ones(p,1)];
        
        f = sparse(1+Nu+q*p,1);
        f(1) = 1;

        prob.c = f;
        prob.a = [[Fu H;Fl Cont] sparse(neq,q*p);...
          sparse(2*(q-1)*Ncheck*NNE,1) I_NM Crit;...
          sparse(p,1+Nu) kron(speye(p),[1 sparse(1,q-1)])];

        prob .blc = beq;
        prob .buc = beq;
        prob .blx = [];
        prob .bux = [];

        prob.cones = cell(p ,1);
        for t=1:p
            prob.cones{t}.type = 'MSK_CT_QUAD';
            prob.cones{t}.sub = 1+Nu+q*(t-1)+(1:q);
        end
           
   
        
        case 'MCt'
        if not(isfield(crit,'red_fact'))
            kred = ones(1,crit.ngz);
        else
            kred = crit.red_fact;
        end
        if not(isfield(crit,'red_fact_t'))
            kred_t = kred;
        else
            kred_t = crit.red_fact_t;
        end
               
        ngz = crit.ngz;
        [az,wz] = quadrature_1D(ngz,'uniform');
        
        q = 9;
        p = ngz*Ncheck*NNE;
        
        
        AN = kron(2/thick.*wz,[speye(3) sparse(3,6)]);
        AM = kron(-wz.*az,[speye(3) sparse(3,6)]);
        Crit = kron(speye(Ncheck*NNE),[AN;AM]);
        
        
        NN = sparse(Ncheck*6,33);
        
        NN(1:6:6*Ncheck,1:3:9) = NNV(1:Ncheck,:);
        NN(2:6:6*Ncheck,2:3:9) = NNV(1:Ncheck,:);
        NN(3:6:6*Ncheck,3:3:9) = NNV(1:Ncheck,:);
        NN(4:6:6*Ncheck,9+(1:3:18)) = NM(1:Ncheck,:);
        NN(5:6:6*Ncheck,9+(2:3:18)) = NM(1:Ncheck,:);
        NN(6:6:6*Ncheck,9+(3:3:18)) = NM(1:Ncheck,:);
        
        I_NM = kron(speye(NNE),NN);        
        

        % rhs members
        beq = sparse(neq+6*Ncheck*NNE+6*p,1);
        beq((neq+6*Ncheck*NNE+1):6:end) = crit.ft.*repmat(kred_t',Ncheck*NNE,1);
        beq((neq+6*Ncheck*NNE+2):6:end) = crit.ft.*repmat(kred_t',Ncheck*NNE,1);
        beq((neq+6*Ncheck*NNE+4):6:end) = crit.fc.*repmat(kred',Ncheck*NNE,1);
        beq((neq+6*Ncheck*NNE+5):6:end) = crit.fc.*repmat(kred',Ncheck*NNE,1);
        
        
        f = sparse(1+Nu+q*p,1);
        f(1) = 1;

        prob.c = f;
        prob.a = [[Fu H;Fl Cont] sparse(neq,q*p);...
          sparse(6*Ncheck*NNE,1) I_NM -Crit;...
          sparse(6*p,1+Nu) kron(speye(p),[[diag([1,1,-sqrt(2)]);diag([-1,-1,-sqrt(2)])] speye(6)])];

        prob .blc = beq;
        prob .buc = beq;
        prob .blx = [];
        prob .bux = [];

        prob.cones = cell(2*p,1);
        for t=1:p
            prob.cones{t}.type = 'MSK_CT_RQUAD';
            prob.cones{t}.sub = 1+Nu+q*(t-1)+(4:6);
            prob.cones{p+t}.type = 'MSK_CT_RQUAD';
            prob.cones{p+t}.sub = 1+Nu+q*(t-1)+(7:9);
        end
        
        
        
        case 'Tresca'
        if not(isfield(crit,'red_fact'))
            kred = ones(1,crit.ngz);
        else
            kred = crit.red_fact;
        end
        if not(isfield(crit,'red_fact_t'))
            kred_t = kred;
        else
            kred_t = crit.red_fact_t;
        end
               
        ngz = crit.ngz;
        [az,wz] = quadrature_1D(ngz,'uniform');
        
        q = 12;
        p = ngz*Ncheck*NNE;
        
        
        AN = kron(2/thick.*wz,[speye(3) sparse(3,q-3)]);
        AM = kron(-wz.*az,[speye(3) sparse(3,q-3)]);
        Crit = kron(speye(Ncheck*NNE),[AN;AM]);
        
        
        NN = sparse(Ncheck*6,33);
        
        NN(1:6:6*Ncheck,1:3:9) = NNV(1:Ncheck,:);
        NN(2:6:6*Ncheck,2:3:9) = NNV(1:Ncheck,:);
        NN(3:6:6*Ncheck,3:3:9) = NNV(1:Ncheck,:);
        NN(4:6:6*Ncheck,9+(1:3:18)) = NM(1:Ncheck,:);
        NN(5:6:6*Ncheck,9+(2:3:18)) = NM(1:Ncheck,:);
        NN(6:6:6*Ncheck,9+(3:3:18)) = NM(1:Ncheck,:);
        
        I_NM = kron(speye(NNE),NN);        
        

        % rhs members
        beq = sparse(neq+6*Ncheck*NNE+9*p,1);
        beq((neq+6*Ncheck*NNE+1):9:end) = repmat(kred_t',Ncheck*NNE,1);
        beq((neq+6*Ncheck*NNE+2):9:end) = repmat(kred_t',Ncheck*NNE,1);
        beq((neq+6*Ncheck*NNE+4):9:end) = repmat(kred',Ncheck*NNE,1);
        beq((neq+6*Ncheck*NNE+5):9:end) = repmat(kred',Ncheck*NNE,1);
        beq((neq+6*Ncheck*NNE+7):9:end) = repmat(kred',Ncheck*NNE,1);
        
        
        f = sparse(1+Nu+q*p,1);
        f(1) = 1;

        prob.c = f;
        prob.a = [[Fu H;Fl Cont] sparse(neq,q*p);...
          sparse(6*Ncheck*NNE,1) I_NM -Crit;...
          sparse(9*p,1+Nu) kron(speye(p),[[diag([1,1,-sqrt(2)]);diag([-1,-1,-sqrt(2)]);[sparse(1,3);0 1 -1;0 0 2]] speye(9)])];

        prob .blc = beq;
        prob .buc = beq;
        prob .blx = [];
        prob .bux = [];

        prob.cones = cell(3*p,1);
        for t=1:p
            prob.cones{t}.type = 'MSK_CT_RQUAD';
            prob.cones{t}.sub = 1+Nu+q*(t-1)+(4:6);
            prob.cones{p+t}.type = 'MSK_CT_RQUAD';
            prob.cones{p+t}.sub = 1+Nu+q*(t-1)+(7:9);
            prob.cones{2*p+t}.type = 'MSK_CT_QUAD';
            prob.cones{2*p+t}.sub = 1+Nu+q*(t-1)+(10:12);
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
        if not(isfield(crit,'red_fact_t'))
            ktred = kcred;
        else
            ktred = crit.red_fact_t;
        end
        
        
        ngz = crit.ngz;
        nga = crit.ac.nbac;
        
        [az,wz] = quadrature_1D(ngz,'uniform');
        xia = crit.ac.zs;
        
        q = 9;
        p = ngz*Ncheck*NNE;
        
        p2 = nga*Ncheck*NNE;
        
        
        AN = kron(2/thick.*wz,[speye(3) sparse(3,6)]);
        AM = kron(-wz.*az,[speye(3) sparse(3,6)]);
        Crit = kron(speye(Ncheck*NNE),[AN;AM]);
        
        Dir = sparse(3,nga);
        for a=1:nga
            Dir(:,a) = [crit.ac.dir(a,1)^2;crit.ac.dir(a,2)^2;crit.ac.dir(a,1)*crit.ac.dir(a,2)];
        end
        AN_ac = repmat(4/thick^2,3,nga).*Dir;
        AM_ac = repmat(-2/thick*xia,3,1).*Dir;
        Crit_ac = kron(speye(Ncheck*NNE),[AN_ac;AM_ac]);
        
        NN = sparse(Ncheck*6,33);
        
        NN(1:6:6*Ncheck,1:3:9) = NNV(1:Ncheck,:);
        NN(2:6:6*Ncheck,2:3:9) = NNV(1:Ncheck,:);
        NN(3:6:6*Ncheck,3:3:9) = NNV(1:Ncheck,:);
        NN(4:6:6*Ncheck,9+(1:3:18)) = NM(1:Ncheck,:);
        NN(5:6:6*Ncheck,9+(2:3:18)) = NM(1:Ncheck,:);
        NN(6:6:6*Ncheck,9+(3:3:18)) = NM(1:Ncheck,:);
        
        I_NM = kron(speye(NNE),NN);        
        

        % rhs members
        beq = sparse(neq+6*Ncheck*NNE+6*p,1);
        beq((neq+6*Ncheck*NNE+1):6:end) = crit.ft.*repmat(ktred',Ncheck*NNE,1);
        beq((neq+6*Ncheck*NNE+2):6:end) = crit.ft.*repmat(ktred',Ncheck*NNE,1);
        beq((neq+6*Ncheck*NNE+4):6:end) = crit.fc.*repmat(kcred',Ncheck*NNE,1);
        beq((neq+6*Ncheck*NNE+5):6:end) = crit.fc.*repmat(kcred',Ncheck*NNE,1);
        
        
        f = sparse(1+Nu+q*p+p2,1);
        f(1) = 1;

        prob.c = f;
        prob.a = [[Fu H;Fl Cont] sparse(neq,q*p+p2);...
          sparse(6*Ncheck*NNE,1) I_NM -Crit -Crit_ac;...
          sparse(6*p,1+Nu) kron(speye(p),[[diag([1,1,-sqrt(2)]);diag([-1,-1,-sqrt(2)])] speye(6)]) sparse(6*p,p2) ];

        prob .blc = beq;
        prob .buc = beq;
        prob .blx = [-inf*ones(1+Nu+q*p,1);-crit.ac.Ns.*repmat(kared',Ncheck*NNE,1)];
        prob .bux = -prob.blx;

        prob.cones = cell(2*p,1);
        for t=1:p
            prob.cones{t}.type = 'MSK_CT_RQUAD';
            prob.cones{t}.sub = 1+Nu+q*(t-1)+(4:6);
            prob.cones{p+t}.type = 'MSK_CT_RQUAD';
            prob.cones{p+t}.sub = 1+Nu+q*(t-1)+(7:9);
        end
        
        
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
        
        
        q = 7;
        C = mat_TC_ell(crit.param_ell,6);
        nell = size(crit.param_ell,2);       
        
        nga = crit.ac.nbac;
        xia = crit.ac.zs;
        
        p = nell*Ncheck*NNE;
        
        p2 = nga*Ncheck*NNE;
        
        
        Crit = kron(speye(Ncheck*NNE),C);
        
        Dir = sparse(3,nga);
        for a=1:nga
            Dir(:,a) = [crit.ac.dir(a,1)^2;crit.ac.dir(a,2)^2;crit.ac.dir(a,1)*crit.ac.dir(a,2)];
        end
        AN_ac = repmat(4/thick^2,3,nga).*Dir;
        AM_ac = repmat(-2/thick*xia,3,1).*Dir;
        Crit_ac = kron(speye(Ncheck*NNE),[AN_ac;AM_ac]);
        
        NN = sparse(Ncheck*6,33);
        
        NN(1:6:6*Ncheck,1:3:9) = NNV(1:Ncheck,:);
        NN(2:6:6*Ncheck,2:3:9) = NNV(1:Ncheck,:);
        NN(3:6:6*Ncheck,3:3:9) = NNV(1:Ncheck,:);
        NN(4:6:6*Ncheck,9+(1:3:18)) = NM(1:Ncheck,:);
        NN(5:6:6*Ncheck,9+(2:3:18)) = NM(1:Ncheck,:);
        NN(6:6:6*Ncheck,9+(3:3:18)) = NM(1:Ncheck,:);
        
        
        I_NM = kron(speye(NNE),NN);        
        

        % rhs members
        beq = sparse(neq+q*Ncheck*NNE,1);
        beq(neq+6*Ncheck*NNE+(1:Ncheck*NNE)) = 1;
        
        
        f = sparse(1+Nu+q*p+p2,1);
        f(1) = 1;

        prob.c = f;
        prob.a = [[Fu H;Fl Cont] sparse(neq,q*p+p2);...
          sparse(6*Ncheck*NNE,1) I_NM -Crit -Crit_ac;...
          sparse(Ncheck*NNE,1+Nu) kron(speye(Ncheck*NNE),repmat([1 sparse(1,q-1)],1,nell)) sparse(Ncheck*NNE,p2)];

        prob .blc = beq;
        prob .buc = beq;
        prob .blx = [-inf*ones(1+Nu+q*p,1);-crit.ac.Ns.*repmat(kared',Ncheck*NNE,1)];
        prob .bux = -prob.blx;

        prob.cones = cell(p,1);
        for t=1:p
            prob.cones{t}.type = 'MSK_CT_QUAD';
            prob.cones{t}.sub = 1+Nu+q*(t-1)+(1:q);
        end
        
       case 'CHell'
      
        
        q = 7;
        C = mat_TC_ell(crit.param_ell,6);
        nell = size(crit.param_ell,2);   
        
        p = nell*Ncheck*NNE;
        
        
        
        Crit = kron(speye(Ncheck*NNE),C);
                
        NN = sparse(Ncheck*6,33);
        
        NN(1:6:6*Ncheck,1:3:9) = NNV(1:Ncheck,:);
        NN(2:6:6*Ncheck,2:3:9) = NNV(1:Ncheck,:);
        NN(3:6:6*Ncheck,3:3:9) = NNV(1:Ncheck,:);
        NN(4:6:6*Ncheck,9+(1:3:18)) = NM(1:Ncheck,:);
        NN(5:6:6*Ncheck,9+(2:3:18)) = NM(1:Ncheck,:);
        NN(6:6:6*Ncheck,9+(3:3:18)) = NM(1:Ncheck,:);
        
        I_NM = kron(speye(NNE),NN);        
        

        % rhs members
        beq = sparse(neq+q*Ncheck*NNE,1);
        beq(neq+6*Ncheck*NNE+(1:Ncheck*NNE)) = 1;
        
        
        f = sparse(1+Nu+q*p,1);
        f(1) = 1;

        prob.c = f;
        prob.a = [[Fu H;Fl Cont] sparse(neq,q*p);...
          sparse(6*Ncheck*NNE,1) I_NM -Crit;
          sparse(Ncheck*NNE,1+Nu) kron(speye(Ncheck*NNE),repmat([1 sparse(1,q-1)],1,nell))];

        prob .blc = beq;
        prob .buc = beq;
        prob.blx = [];
        prob.bux = [];

        prob.cones = cell(p,1);
        for t=1:p
            prob.cones{t}.type = 'MSK_CT_QUAD';
            prob.cones{t}.sub = 1+Nu+q*(t-1)+(1:q);
        end
        
    case 'Ilyushin'
        if not(isfield(crit,'red_fact'))
            kred = ones(1,crit.ngz);
        else
            kred = crit.red_fact;
        end
        
                       
        
        
        q = 7;
        p = 2*Ncheck*NNE;
        
        P = sparse([1 -1/2 0;0 sqrt(3)/2 0; 0 0 sqrt(3)]);
        P1 = [P sqrt(3)/6.*P; sparse(3,3) sqrt(33)/6.*P]*diag([thick/4*ones(1,3),ones(1,3)]);
        P2 = [P -sqrt(3)/6.*P; sparse(3,3) sqrt(33)/6.*P]*diag([thick/4*ones(1,3),ones(1,3)]);
        
        Crit = kron(speye(Ncheck*NNE),[P1;P2]);
        
        
        NN = sparse(Ncheck*6,33);
        
        NN(1:6:6*Ncheck,1:3:9) = NNV(1:Ncheck,:);
        NN(2:6:6*Ncheck,2:3:9) = NNV(1:Ncheck,:);
        NN(3:6:6*Ncheck,3:3:9) = NNV(1:Ncheck,:);
        NN(4:6:6*Ncheck,9+(1:3:18)) = NM(1:Ncheck,:);
        NN(5:6:6*Ncheck,9+(2:3:18)) = NM(1:Ncheck,:);
        NN(6:6:6*Ncheck,9+(3:3:18)) = NM(1:Ncheck,:);
        
        I_NM = kron(speye(NNE),NN);        
        

        % rhs members
        beq = [sparse(neq+(q-1)*p,1);ones(p,1)];
        
        f = sparse(1+Nu+q*p,1);
        f(1) = 1;

        prob.c = f;
        prob.a = [[Fu H;Fl Cont] sparse(neq,q*p);...
          sparse((q-1)*p,1) Crit*I_NM kron(speye(p),[sparse(q-1,1) speye(q-1)]);...
          sparse(p,1+Nu) kron(speye(p),[1 sparse(1,q-1)])];

        prob .blc = beq;
        prob .buc = beq;
        prob .blx = [];
        prob .bux = [];

        prob.cones = cell(p ,1);
        for t=1:p
            prob.cones{t}.type = 'MSK_CT_QUAD';
            prob.cones{t}.sub = 1+Nu+q*(t-1)+(1:q);
        end
           
        
    otherwise
        error('TODO')
end

end

