function Cont = continuity(mesh)
%CONTINUITY Summary of this function goes here
%   Detailed explanation goes here
NED_rt = mesh.act_edges.rottan.NED;
NED_rn = mesh.act_edges.rotnorm.NED;
NED_utr = mesh.act_edges.utrans.NED;
NED_utg = mesh.act_edges.utang.NED;
NED_un = mesh.act_edges.unorm.NED;

NEQ = mesh.NEQ;
Cont_rt = sparse(3*NED_rt,NEQ);
Cont_rn = sparse(3*NED_rn,NEQ);
Cont_utr = sparse(2*NED_utr,NEQ);
Cont_utg = sparse(2*NED_utg,NEQ);
Cont_un = sparse(2*NED_un,NEQ);

sign_idrn = ones(NED_rn,1);
sign_idrt = ones(NED_rt,1);
sign_idutr = ones(NED_utr,1);
sign_idutg = ones(NED_utg,1);
sign_idun = ones(NED_un,1);

for e=1:mesh.NNE
    node = mesh.connec(e,:);
    T = mesh.coor(node,:);
    
    i1 = [1 2 3];
    i2 = [2 3 1];
    i3 = [4 5 6];
    for j=1:3
        tang = T(i2(j),:)'-T(i1(j),:)';
        tang = tang./norm(tang);
        
        Rloc = rep_loc(T);
        a1 = Rloc(:,1);
        a2 = Rloc(:,2);
        zloc = Rloc(:,3);
        
        if (node(i1(j))>node(i2(j)));
            tang = -tang;
            temp = i1(j);
            i1(j) = i2(j);
            i2(j) = temp;
        end
        
        id_bc = mesh.bc_edges.node2edg(node(i1(j)),node(i2(j)));
        
        id = mesh.act_edges.all.node2edg(node(i1(j)),node(i2(j)));
        zloc_av = mesh.edg_norm(id,:)';
%         a1 = cross(a2,zloc);
%         a2 = cross(zloc,a1);
        
        if (id_bc >0) % if edge is on the boundary, the normal must point in the outward direction
            xG = mean(T(1:3,:));
            xM = mean(T([i1(j) i2(j)],:));
            normal = cross(tang,zloc);
            normal = sign((xM-xG)*normal).*normal;
            normal_av = cross(tang,zloc_av);
            normal_av = sign((xM-xG)*normal_av).*normal_av;
            s = 1;
        else
            xG = mean(T(1:3,:));
            xM = mean(T([i1(j) i2(j)],:));
            normal = cross(tang,zloc);
            normal_av = cross(tang,zloc_av);
            % s = 1;
            s = sign((xM-xG)*normal_av);
        end
        
        % expression of M_{nn},M_{nt} from their components M_{ab} in the
        % local basis of the element
        Mnn = [(normal'*a1)*(normal_av'*a1) (normal'*a2)*(normal_av'*a2) (normal'*a1)*(normal_av'*a2)+(normal_av'*a1)*(normal'*a2)];
        Mnt = [(normal'*a1)*(tang'*a1) (normal'*a2)*(tang'*a2) (normal'*a1)*(tang'*a2)+(normal'*a2)*(tang'*a1)];
        
        ProjM = sparse(6,18);
        ProjM(1,3*(i1(j)-1)+(1:3)) = Mnn;
        ProjM(2,3*(i2(j)-1)+(1:3)) = Mnn;
        ProjM(3,3*(i3(j)-1)+(1:3)) = Mnn;
        ProjM(4,3*(i1(j)-1)+(1:3)) = Mnt;
        ProjM(5,3*(i2(j)-1)+(1:3)) = Mnt;
        ProjM(6,3*(i3(j)-1)+(1:3)) = Mnt;
        
        nu1 = normal'*a1;
        nu1_av = normal_av'*a1;
        t1 = tang'*a1;
        z1 = zloc'*a1;
        z1_av = zloc_av'*a1;
        nu2 = normal'*a2;
        nu2_av = normal_av'*a2;
        t2 = tang'*a2;
        z2 = zloc'*a2;
        z2_av = zloc_av'*a2;
        ProjR = sparse(6,15);
        l1 = [3*(i1(j)-1)+(1:3) 9+2*(i1(j)-1)+(1:2)];
        ProjR(1,l1) = [nu1*nu1_av nu2*nu2_av nu2*nu1_av+nu2_av*nu1 (zloc'*normal_av)*nu1 (zloc'*normal_av)*nu2];
        ProjR(2,l1) = [nu1*t1 nu2*t2 nu2*t1+nu1*t2 (zloc'*tang)*nu1 (zloc'*tang)*nu2];
        ProjR(3,l1) = [nu1*z1_av nu2*z2_av nu2*z1_av+nu1*z2_av nu1*(zloc'*zloc_av) nu2*(zloc'*zloc_av)];
        
        l2 = [3*(i2(j)-1)+(1:3) 9+2*(i2(j)-1)+(1:2)];
        ProjR(4:6,l2) = ProjR(1:3,l1);
        
        id_rn = mesh.act_edges.rotnorm.node2edg(node(i1(j)),node(i2(j)));
        id_rt = mesh.act_edges.rottan.node2edg(node(i1(j)),node(i2(j)));
        id_utr = mesh.act_edges.utrans.node2edg(node(i1(j)),node(i2(j)));
        id_utg = mesh.act_edges.utang.node2edg(node(i1(j)),node(i2(j)));
        id_un = mesh.act_edges.unorm.node2edg(node(i1(j)),node(i2(j)));
        
        if (s<0)
            sign_idrn(id_rn) = -1;
            sign_idrt(id_rt) = -1;
            sign_idutr(id_utr) = -1;
            sign_idutg(id_utg) = -1;
            sign_idun(id_un) = -1;
        end
        
        leM = 33*(e-1)+9+(1:18);
        leNV = 33*(e-1)+ [1:9 28:33];
        
        if id_rn>0
            Cont_rn(3*(id_rn-1)+(1:3),leM) = Cont_rn(3*(id_rn-1)+(1:3),leM) + sign_idrn(id_rn).*ProjM(1:3,:); 
            sign_idrn(id_rn)=-sign_idrn(id_rn);
        end
        if id_rt>0
            Cont_rt(3*(id_rt-1)+(1:3),leM) = Cont_rt(3*(id_rt-1)+(1:3),leM) + sign_idrt(id_rt).*ProjM(4:6,:);
            sign_idrt(id_rt)=-sign_idrt(id_rt);
        end
        if id_utr>0
            Cont_utr(2*(id_utr-1)+(1:2),leNV) = Cont_utr(2*(id_utr-1)+(1:2),leNV) + sign_idutr(id_utr).*ProjR([3 6],:);
            sign_idutr(id_utr)=-sign_idutr(id_utr);
        end
        if id_utg>0
            Cont_utg(2*(id_utg-1)+(1:2),leNV) = Cont_utg(2*(id_utg-1)+(1:2),leNV) + sign_idutg(id_utg).*ProjR([2 5],:);
            sign_idutg(id_utg)=-sign_idutg(id_utg);
        end
        if id_un>0
            Cont_un(2*(id_un-1)+(1:2),leNV) = Cont_un(2*(id_un-1)+(1:2),leNV) + sign_idun(id_un).*ProjR([1 4],:);
            sign_idun(id_un)=-sign_idun(id_un);
        end
    end
    
end

Cont = [Cont_rn;Cont_rt;Cont_un;Cont_utg;Cont_utr;];
end

