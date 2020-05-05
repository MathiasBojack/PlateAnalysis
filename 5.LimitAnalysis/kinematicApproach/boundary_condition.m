function BC = boundary_condition(mesh)
%boundary_condition Summary of this function goes here
%   normal and tangent are related to element edge (not element normal)
% Boundary conditions are expressed on normal inplane displacement,
%   tangential inplane displacement, transversal displacement,
%   rotation along tangent, rotation around normal, rotation around z
global dim

NEQ = mesh.NEQ;
BC = [];
buff = 0;
for i_ed=1:mesh.nb_edg
    bc_id = logical([mesh.bc.utang(i_ed);mesh.bc.unorm(i_ed);mesh.bc.utrans(i_ed);0;mesh.bc.utrans(i_ed);1]);
    n_bc = sum(bc_id);
    
    if (n_bc>0)
        pos_edg = find(mesh.edges(:,3)==i_ed);
        %nd_edg = unique(mesh.nd_edg_list{i_ed});
        nd_edg = unique(mesh.edges(pos_edg,1:2));
        BC_tmp = sparse(n_bc*length(nd_edg),NEQ);
        for e = 1:length(pos_edg)
            elem = find_el_from_node(mesh.edges(pos_edg,1:2),mesh.connec);
            nodes = mesh.connec(elem(e),:);
            i1 = find(nodes == mesh.edges(pos_edg(e),1));
            i2 = find(nodes == mesh.edges(pos_edg(e),2));
            T = mesh.coor(mesh.connec(elem(e),:),:);
            Rloc = rep_loc(T);
            tang = T(i2,:)' - T(i1,:)';
            tang = tang./norm(tang);
            normal = cross(Rloc(:,3),tang);
            Rtn = [tang normal Rloc(:,3)];

            p1 = n_bc*(find(nd_edg==nodes(i1))-1) + (1:n_bc);
            p2 = n_bc*(find(nd_edg==nodes(i2))-1) + (1:n_bc);

            Rtot = [Rtn' zeros(3,3);zeros(3,3) Rtn'];

            pe1 = dim*(nodes(i1)-1) + (1:dim);
            pe2 = dim*(nodes(i2)-1) + (1:dim);

            BC_tmp(p1,pe1) = BC_tmp(p1,pe1) + Rtot(bc_id,:);
            BC_tmp(p2,pe2) = BC_tmp(p2,pe2) + Rtot(bc_id,:);
        end
    end
    
    BC= [BC;BC_tmp];
end

% IP_rot = sparse(3*mesh.NNO,NEQ);
% for e=1:mesh.NNE
%     node = mesh.connec(e,:);
%     T = mesh.coor(node,:);
%     Rloc = rep_loc(T);
%     Rot = rot_loc2glob(Rloc);
%     pe = LtoGmap(node,dim*mesh.NNO+e);
%     idx = [3;9;15];
%     %idx = [1 3;7 9; 13 15];
%     n = size(idx,2);
%     for i=1:3
%         IP_rot(n*(node(i)-1)+(1:n),pe) = IP_rot(n*(node(i)-1)+(1:n),pe) + Rot(:,idx(i,:))';
%     end
% end
% BC= [BC;IP_rot];


L = sparse(3,NEQ);
L(1,1) = 1;
L(2,2) = 1;
L(3,3) = 1;
% 
% L(2,dim*(1-1) +2) = 1;
% L(2,dim*(3-1) +2)= -1;
BC = [BC;L];
end
    
%     
%     if (bcu(e)==1)                      % for boundaries with zero displacement
%     pos_e = mesh.edges(:,3) == e;
%     AA = mesh.edges(pos_e,1:2);
%     A = [mesh.edges(pos_e,1);mesh.edges(pos_e,2)];
%     nb_e = size(AA,1);
%     
%     [b nn mm] = unique(A);
%     for k=1:size(b,1)
%        curr_node = b(k);
%        pk = find(mm == k);
%        ppk = mod(pk-1,nb_e)+1;
%        neighb_edg = AA(ppk,:);
%        
%        if (size(neighb_edg,1)==1)
%            T = coor(neighb_edg,:);
%            L = sqrt((T(2,1)-T(1,1))^2+(T(2,2)-T(1,2))^2);
%            tang = (T(2,:)-T(1,:))/L;               % tangent vector of the current edge
%       
%        elseif (size(neighb_edg,1)==2)
%            T1 = coor(neighb_edg(1,:)',:);
%            L1 = sqrt((T1(2,1)-T1(1,1))^2+(T1(2,2)-T1(1,2))^2);
%            tang1 = (T1(2,:)-T1(1,:))/L1;               % tangent vector of the first neighbouring edge
%            T2 = coor(neighb_edg(2,:)',:);
%            L2 = sqrt((T2(2,1)-T2(1,1))^2+(T2(2,2)-T2(1,2))^2);
%            tang2 = (T2(2,:)-T2(1,:))/L2;               % tangent vector of the second neighbouring edge
%            
%            tang = (tang1+tang2)./2;                    % averaging tangent vector
%            tang = tang./norm(tang);
%        else
%            error('Problem in the neighbouring edges calculation')
%        end
%            buff = buff + 1;
%            BC(buff,dof(dim*(curr_node-1)+(2:3))) = BC(buff,dof(dim*(curr_node-1)+(2:3))) + tang;
%    
%     end
%     end
% end
% 
% 
% BC = BC(1:buff,:);
% end

