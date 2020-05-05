function bc_type = active_edges(mesh,bc)
%ACTIVE_EDGES : tables for edges which are active for the given boundary
% condition (bc)
% bc_type.NED = number of active edges
% bc_type.node2edg(n1,n2) = id   : matrix giving the global numbering (id)
% of the edge having n1 and n2 for end nodes
% bc_type.edg2node(id) = [n1 n2] : reciprocal of node2edg

NNO = mesh.NNO;
NNE = mesh.NNE;

node2edg = sparse(NNO,NNO);
edge_bc = find(bc==1);
buffn=0;
for i = 1:length(edge_bc)
   edg = edge_bc(i);
   pos_nd_edg = mesh.edges(:,3)==edg;
   nd_edg = mesh.edges(pos_nd_edg,1:2);
   if (size(nd_edg)~=[0 0])
       for ii=1:size(nd_edg,1)
       buffn=buffn-1;
       node2edg(nd_edg(ii,1),nd_edg(ii,2))=buffn;
       node2edg(nd_edg(ii,2),nd_edg(ii,1))=buffn;
       end
   end
end

buffp=0;
edg2node=[];
for e=1:NNE
   connec = mesh.connec(e,:);
   nnode = 3;
   i1 = 1:nnode;
   i2 = [2:nnode 1];
   for i=1:nnode
       ii1 = connec(i1(i));
       ii2 = connec(i2(i));
       if (node2edg(ii1,ii2)==0)
            buffp = buffp+1;
            node2edg(ii1,ii2)=buffp;
            node2edg(ii2,ii1)=buffp;
            edg2node = [edg2node; ii1 ii2];
       end
   end
end

bc_type.NED = buffp;
bc_type.node2edg = node2edg;
bc_type.edg2node = edg2node;
end

