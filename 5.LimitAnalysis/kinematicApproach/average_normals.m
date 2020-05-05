function mesh = average_normals(mesh)
%AVERAGE_NORMALS Summary of this function goes here
%   Detailed explanation goes here

edg_norm = zeros(mesh.NED,3);
edg = zeros(mesh.NED,3);
node_norm = zeros(mesh.NNO,3);

for e=1:mesh.NNE
    node = mesh.connec(e,:);
    T = mesh.coor(node,:);
    Rloc = rep_loc(T);
    
    i1 = [1 2 3];
    i2 = [2 3 1];
    for j=1:3
         id = mesh.act_edges.prm.node2edg(node(i1(j)),node(i2(j)));
         node_norm(node(i1(j)),:) = node_norm(node(i1(j)),:) + Rloc(:,3)';
         node_norm(node(i2(j)),:) = node_norm(node(i2(j)),:) + Rloc(:,3)';
         
         if (id>0)
%            edg_norm(id,:) = edg_norm(id,:) + Rloc(:,3)';
         
            edg(id,:) = mesh.coor(node(i2(j)),:) - mesh.coor(node(i1(j)),:);
            edg(id,:) = edg(id,:)./norm(edg(id,:));
            
            elem_at_disc = mesh.id_dof.elem_at_disc{id};
            if length(elem_at_disc)==1
                pe1 = 1;
            else
                pe1 = find(elem_at_disc(1:end-1)==e);
            end
            
            if not(isempty(pe1))
                id_g1 = sum(mesh.id_dof.nb_disc(1:id-1))+pe1;     
                edg_norm(id_g1,:) = edg_norm(id_g1,:)+Rloc(:,3)';
            end

            if length(elem_at_disc)>1
                pe2 = find(elem_at_disc(2:end)==e);
                if not(isempty(pe2))
                    id_g2 = sum(mesh.id_dof.nb_disc(1:id-1))+pe2;
                    edg_norm(id_g2,:) = edg_norm(id_g2,:)+Rloc(:,3)';
                    
                end
            end
         end
         
    end
end
mesh.edg_norm = edg_norm./(sqrt(edg_norm(:,1).^2+edg_norm(:,2).^2+edg_norm(:,3).^2)*[1 1 1]);
mesh.node_norm = node_norm./(sqrt(node_norm(:,1).^2+node_norm(:,2).^2+node_norm(:,3).^2)*[1 1 1]);
mesh.inplane_norm = cross(mesh.edg_norm',edg')';
end

