function [Prm_strain,Prm_disc,diss,diss_renf] = prm_post_process(d,dd,BC,crit,thick,edg_lgth,mesh,ngz,nged,normal_loc)
%BUILDING_MOSEK Construction of the Mosek optimization problems

mesh.thick_list  = zeros(mesh.NNE,1);
for j=1:length(thick)
    mesh.thick_list(mesh.el_tag==j) = thick(j);
end
mesh.edg_thick_list = zeros(mesh.NED,1);
for j=1:size(mesh.id_dof.nb_disc,1)
    el = mesh.id_dof.elem_at_disc{j};
    loc_thick = mesh.thick_list(el);
    jj = sum(mesh.id_dof.nb_disc(1:j-1)) +(1:mesh.id_dof.nb_disc(j)) ;
    if length(jj)>1
        mesh.edg_thick_list(jj) = min(loc_thick(1:end-1),loc_thick(2:end));
    else
        mesh.edg_thick_list(jj) = min(loc_thick);
    end
end

NNE = mesh.NNE;
NED = mesh.NED;
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
[az,wz] = quadrature_1D(ngz,'trapeze');

switch crit.type
                
    case 'vonMises'
        if not(isfield(crit,'red_fact'))
            kred = ones(1,ngz);
        else
            kred = crit.red_fact;
        end
        


        
        q1 = 4; % number of conic variables for Prm_strain/point
        q2 = 3; % number of conic variables for Prm_disc/point
        p1 = ngz*3*NNE; % total number of integration points for Prm_strain        
        Nd = mesh.Nu;
        
        
        Cb = sparse([2 1 0;0 sqrt(3) 0;0 0 1]);
        Cn = sparse(diag([2,1]));
        
        p2 = ngz*nged*NED; % total number of intergration points for Prm_disc

        c_strain = repmat(kron(kred'.*wz',(crit.strength.*area)'./3),3,1);
        c_disc = kron(wd',kron(kred'.*wz',crit.strength_ed'.*edg_lgth'./2));
        
        
        N_edlin = [(1-ad')/2 (ad'+1)./2];
        N_edquad = [(ad'-1).*ad'/2 (ad'+1).*ad'/2 1-ad'.^2];
                
        Cn_r = kron(-az',thick^2/4*Cn./sqrt(3));
        Cn_u = kron(ones(ngz,1),[1/2*thick*Cn./sqrt(3),sparse(2,1)]);
        Crit_disc = [kron(N_edquad,Cn_u) ...
                        kron(N_edlin,Cn_r) sparse(ngz*2*nged,2)];
        
        Crit_c = kron(-az',thick^2/4*Cb./sqrt(3));
        Crit_e = kron(ones(ngz,1),1/2*thick*Cb./sqrt(3));
                    
        Crit_strain = [repmat(Crit_c,3,1) kron(speye(3),Crit_e) sparse(9*ngz,6)];
        C_disc = kron(speye(NED),Crit_disc);

        ll = setdiff((1:sd*NED)',BC);
        r_strain = kron(speye(NNE),Crit_strain)*d;
        r_disc = C_disc(:,ll)*dd;
        
        t_strain = sqrt(sum(reshape(r_strain.^2,q1-1,p1))');
        t_disc = sqrt(sum(reshape(r_disc.^2,q2-1,p2))');
        
        
        Prm_strain = c_strain(:)'*t_strain;
        Prm_disc = c_disc(:)'*t_disc;
            
                
    case 'vonMises-varthick'
        if not(isfield(crit,'red_fact'))
            kred = ones(1,ngz);
        else
            kred = crit.red_fact;
        end
        


        
        q1 = 4; % number of conic variables for Prm_strain/point
        q2 = 3; % number of conic variables for Prm_disc/point
        p1 = ngz*3*NNE; % total number of integration points for Prm_strain        
        Nd = mesh.Nu;
        
        
        Cb = sparse([2 1 0;0 sqrt(3) 0;0 0 1]);
        Cn = sparse(diag([2,1]));
        
        p2 = ngz*nged*NED; % total number of intergration points for Prm_disc

        c_strain = repmat(kron(kred'.*wz',(crit.strength.*area)'./3),3,1);
        c_disc = kron(wd',kron(kred'.*wz',crit.strength_ed'.*edg_lgth'./2));
        
        
        N_edlin = [(1-ad')/2 (ad'+1)./2];
        N_edquad = [(ad'-1).*ad'/2 (ad'+1).*ad'/2 1-ad'.^2];
                
        Cn_r = kron(-az',Cn./sqrt(3));
        Cn_u = kron(ones(ngz,1),[1/2*Cn./sqrt(3),sparse(2,1)]);
        Crit_disc = [kron(N_edquad,Cn_u) ...
                        kron(N_edlin,Cn_r) sparse(ngz*2*nged,2)];
        
        Crit_c = kron(-az',Cb./sqrt(3));
        Crit_e = kron(ones(ngz,1),1/2*Cb./sqrt(3));
                    
        Ithick_s = kron(spdiags(mesh.thick_list,0,NNE,NNE),[sparse(9*ngz,3) kron(speye(3),kron(ones(ngz,1),ones(3,3))) sparse(9*ngz,6)]);
        Ithick_s = Ithick_s + kron(spdiags(mesh.thick_list.^2/4,0,NNE,NNE),[repmat(kron(ones(ngz,1),speye(3)),3,1) sparse(9*ngz,15)]);
        
        Ithick_d = kron(spdiags(mesh.edg_thick_list,0,NED,NED),[kron(speye(3),kron(ones(ngz,1),ones(2,3))) sparse(6*ngz,6)]);
        Ithick_d = Ithick_d + kron(spdiags(mesh.edg_thick_list.^2/4,0,NED,NED),[sparse(6*ngz,9) kron(ones(3,2),kron(ones(ngz,1),speye(2))) sparse(6*ngz,2)]);
        
        Crit_strain = [repmat(Crit_c,3,1) kron(speye(3),Crit_e) sparse(9*ngz,6)];
        C_disc = kron(speye(NED),Crit_disc);

        ll = setdiff((1:sd*NED)',BC);
        r_strain = (kron(speye(NNE),Crit_strain).*Ithick_s)*d;
        r_disc = C_disc(:,ll)*dd;
        
        t_strain = sqrt(sum(reshape(r_strain.^2,q1-1,p1))');
        t_disc = sqrt(sum(reshape(r_disc.^2,q2-1,p2))');
        
        
        Prm_strain = c_strain(:)'*t_strain;
        Prm_disc = c_disc(:)'*t_disc;
            
    case 'MCt'
        
        if not(isfield(crit,'red_fact'))
            kred = ones(1,ngz);
        else
            kred = crit.red_fact;
        end
        
        

        c_disc = kron(wd',kron(kred'.*wz',edg_lgth'./2));
                
        N_edlin = [(1-ad')/2 (ad'+1)./2];
        N_edquad = [(ad'-1).*ad'/2 (ad'+1).*ad'/2 1-ad'.^2];
                
        Cn = sparse([1 0;0 0;0 1]); 
        Cn_r = kron(-thick^2/4.*az',Cn);
        Cn_u = kron(ones(ngz,1),[thick/2*Cn,sparse(3,1)]);
        Crit_disc = [kron(N_edquad,Cn_u) ...
                        kron(N_edlin,Cn_r) sparse(ngz*3*nged,2)];
        
                    
        Cb = speye(3);           
        Crit_c = kron(-az',thick^2/4.*Cb);
        Crit_e = kron(ones(ngz,1),thick/2*Cb);
                    
        Crit_strain = [repmat(Crit_c,3,1) kron(speye(3),Crit_e) sparse(9*ngz,6)];
        
        C_disc = kron(speye(NED),Crit_disc);

        ll = setdiff((1:sd*NED)',BC);
        r_strain = kron(speye(NNE),Crit_strain)*d;
        r_disc = C_disc(:,ll)*dd;
        
        rs_I = (r_strain(1:3:end)+r_strain(2:3:end))./2 + sqrt((r_strain(1:3:end)-r_strain(2:3:end)).^2+r_strain(3:3:end).^2)./2;
        rs_II = (r_strain(1:3:end)+r_strain(2:3:end))./2 - sqrt((r_strain(1:3:end)-r_strain(2:3:end)).^2+r_strain(3:3:end).^2)./2;
        rd_I = (r_disc(1:3:end)+r_disc(2:3:end))./2 + sqrt((r_disc(1:3:end)-r_disc(2:3:end)).^2+r_disc(3:3:end).^2)./2;
        rd_II = (r_disc(1:3:end)+r_disc(2:3:end))./2 - sqrt((r_disc(1:3:end)-r_disc(2:3:end)).^2+r_disc(3:3:end).^2)./2;
                
        cs_t = repmat(kron(kred'.*wz',(crit.ft.*area)'./3),3,1);
        cs_c = repmat(kron(kred'.*wz',(crit.fc.*area)'./3),3,1);
        cd_t = crit.ft.*c_disc;
        cd_c = crit.fc.*c_disc;
        
        Prm_strain = sum(max(cs_t(:).*rs_I,-cs_c(:).*rs_I)) + sum(max(cs_t(:).*rs_II,-cs_c(:).*rs_II));
        Prm_disc = sum(max(cd_t(:).*rd_I,-cd_c(:).*rd_I)) + sum(max(cd_t(:).*rd_II,-cd_c(:).*rd_II));
        

    case 'Tresca'
        
        if not(isfield(crit,'red_fact'))
            kred = ones(1,ngz);
        else
            kred = crit.red_fact;
        end
        
        
        

        c_disc = kron(wd',kron(kred'.*wz',edg_lgth'./2));
                
        N_edlin = [(1-ad')/2 (ad'+1)./2];
        N_edquad = [(ad'-1).*ad'/2 (ad'+1).*ad'/2 1-ad'.^2];
                
        Cn = sparse([1 0;0 0;0 1]); 
        Cn_r = kron(-thick^2/4.*az',Cn);
        Cn_u = kron(ones(ngz,1),[thick/2*Cn,sparse(3,1)]);
        Crit_disc = [kron(N_edquad,Cn_u) ...
                        kron(N_edlin,Cn_r) sparse(ngz*3*nged,2)];
        
                    
        Cb = speye(3);           
        Crit_c = kron(-az',thick^2/4.*Cb);
        Crit_e = kron(ones(ngz,1),thick/2*Cb);
                    
        Crit_strain = [repmat(Crit_c,3,1) kron(speye(3),Crit_e) sparse(9*ngz,6)];
        
        C_disc = kron(speye(NED),Crit_disc);

        ll = setdiff((1:sd*NED)',BC);
        r_strain = kron(speye(NNE),Crit_strain)*d;
        r_disc = C_disc(:,ll)*dd;
        
        rs_I = (r_strain(1:3:end)+r_strain(2:3:end))./2 + sqrt((r_strain(1:3:end)-r_strain(2:3:end)).^2+r_strain(3:3:end).^2)./2;
        rs_II = (r_strain(1:3:end)+r_strain(2:3:end))./2 - sqrt((r_strain(1:3:end)-r_strain(2:3:end)).^2+r_strain(3:3:end).^2)./2;
        rd_I = (r_disc(1:3:end)+r_disc(2:3:end))./2 + sqrt((r_disc(1:3:end)-r_disc(2:3:end)).^2+r_disc(3:3:end).^2)./2;
        rd_II = (r_disc(1:3:end)+r_disc(2:3:end))./2 - sqrt((r_disc(1:3:end)-r_disc(2:3:end)).^2+r_disc(3:3:end).^2)./2;
                
        cs = repmat(kron(kred'.*wz',(area)'./3),3,1);
        cd = c_disc;
        
        Prm_strain = sum(max(max(abs(cs(:).*rs_I),abs(cs(:).*rs_II)),abs(cs(:).*(rs_I+rs_II))));
        Prm_disc = sum(max(max(abs(cd(:).*rd_I),abs(cd(:).*rd_II)),abs(cd(:).*(rd_I+rd_II))));     
            
    case {'MCt-r','MCt-r2'}
        if not(isfield(crit,'red_fact'))
            kcred = ones(1,ngz);
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
        
        
        nga = crit.ac.nbac;
        
        [ad,wd] = quadrature_1D(nged,'trapeze');
        Nedg_lin = [(1-ad')./2 (1+ad')./2];
        Nedg_quad = [ad'.*(ad'-1)./2 ad'.*(1+ad')./2  1-ad'.^2];
        
        [az,wz] = quadrature_1D(ngz,'trapeze');
        xia = crit.ac.zs;

        c_disc_c = kron(wd',kron(kcred'.*wz',edg_lgth'./2));
        c_disc_t = kron(wd',kron(ktred'.*wz',edg_lgth'./2));
                
        N_edlin = [(1-ad')/2 (ad'+1)./2];
        N_edquad = [(ad'-1).*ad'/2 (ad'+1).*ad'/2 1-ad'.^2];
                
        Cn = sparse([1 0;0 0;0 1]); 
        Cn_r = kron(-thick^2/4*az',Cn);
        Cn_u = kron(ones(ngz,1),[thick/2*Cn,sparse(3,1)]);
        Crit_disc = [kron(N_edquad,Cn_u) ...
                        kron(N_edlin,Cn_r) sparse(ngz*3*nged,2)];
        
                    
        Cb = speye(3);           
        Crit_c = kron(-thick^2/4*az',Cb);
        Crit_e = kron(ones(ngz,1),thick/2*Cb);
                    
        Crit_strain = [repmat(Crit_c,3,1) kron(speye(3),Crit_e) sparse(9*ngz,6)];
        
        C_disc = kron(speye(NED),Crit_disc);

        ll = setdiff((1:sd*NED)',BC);
        r_strain = kron(speye(NNE),Crit_strain)*d;
        r_disc = C_disc(:,ll)*dd;
        
        rs_I = (r_strain(1:3:end)+r_strain(2:3:end))./2 + sqrt((r_strain(1:3:end)-r_strain(2:3:end)).^2+r_strain(3:3:end).^2)./2;
        rs_II = (r_strain(1:3:end)+r_strain(2:3:end))./2 - sqrt((r_strain(1:3:end)-r_strain(2:3:end)).^2+r_strain(3:3:end).^2)./2;
        rd_I = (r_disc(1:3:end)+r_disc(2:3:end))./2 + sqrt((r_disc(1:3:end)-r_disc(2:3:end)).^2+r_disc(3:3:end).^2)./2;
        rd_II = (r_disc(1:3:end)+r_disc(2:3:end))./2 - sqrt((r_disc(1:3:end)-r_disc(2:3:end)).^2+r_disc(3:3:end).^2)./2;
                
        cs_t = repmat(kron(ktred'.*wz',(crit.ft.*area)'./3),3,1);
        cs_c = repmat(kron(kcred'.*wz',(crit.fc.*area)'./3),3,1);
        cd_t = crit.ft.*c_disc_t;
        cd_c = crit.fc.*c_disc_c;
        
        Prm_strain = sum(max(cs_t(:).*rs_I,-cs_c(:).*rs_I)) + sum(max(cs_t(:).*rs_II,-cs_c(:).*rs_II))
        Prm_disc = sum(max(cd_t(:).*rd_I,-cd_c(:).*rd_I)) + sum(max(cd_t(:).*rd_II,-cd_c(:).*rd_II))
        
        
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
        r_renf_strain = kron(speye(NNE),Crit_renf_strain)*d;
        cs_renf = crit.ac.Ns.*repmat(kron(kared',area'./3),3,1);
        Prm_renf_strain = cs_renf(:)'*abs(r_renf_strain)
        
        
        Cn = sparse([1 0]); 
        Cn_u = kron(ones(nga,1),[1 0 0]);
        Crit_disc_n = [kron(sparse(Nedg_quad),Cn_u) ...
                        kron(sparse(Nedg_lin),kron(-xia'*thick/2,Cn)) sparse(nged*nga,2)];
        Ct = sparse([0 1]); 
        Ct_u = kron(ones(nga,1),[0 1 0]);
        Crit_disc_t = [kron(sparse(Nedg_quad),Ct_u) ...
                        kron(sparse(Nedg_lin),kron(-xia'*thick/2,Ct))  sparse(nged*nga,2)];
        C_disc_renf = sparse(diag(Dir_n(:).^2))*kron(speye(NED),Crit_disc_n) + sparse(diag(Dir_n(:).*Dir_t(:)))*kron(speye(NED),Crit_disc_t);
        r_disc_renf = C_disc_renf(:,ll)*dd;
        cd_renf = crit.ac.Ns.*kron(wd',kron(kared',edg_lgth'./2));
        Prm_renf_disc = cd_renf(:)'*abs(r_disc_renf)
        
        diss = max(cs_t(:).*rs_I,-cs_c(:).*rs_I) + max(cs_t(:).*rs_II,-cs_c(:).*rs_II);
        diss = sum(reshape(diss,ngz,3*mesh.NNE))';
        
        diss_renf = sum(reshape(cs_renf(:).*abs(r_renf_strain),nga,3*mesh.NNE))';
        
        Prm_strain = Prm_strain + Prm_renf_strain; 
        Prm_disc = Prm_disc + Prm_renf_disc; 
end  


end

