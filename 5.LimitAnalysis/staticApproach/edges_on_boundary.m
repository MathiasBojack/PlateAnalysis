function bc_edges = edges_on_boundary(mesh)
%ACTIVE_EDGES Summary of this function goes here
%   Detailed explanation goes here
NNO = mesh.NNO;
NNE = mesh.NNE;

node2edg = sparse(NNO,NNO);
buffp=0;

for edg=1:length(mesh.bc.utrans)
   pos_nd_edg = mesh.edges(:,3)==edg;
   nd_edg = mesh.edges(pos_nd_edg,1:2);
    for ii=1:size(nd_edg,1)
       buffp=buffp+1;
       node2edg(nd_edg(ii,1),nd_edg(ii,2))=buffp;
       node2edg(nd_edg(ii,2),nd_edg(ii,1))=buffp;
    end
end

bc_edges.NED = buffp;
bc_edges.node2edg = node2edg;
end

