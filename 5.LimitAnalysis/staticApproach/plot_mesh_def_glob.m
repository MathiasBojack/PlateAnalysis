function [ output_args ] = plot_mesh_def_glob(mesh,U,alp_def)
%PLOT_MESH Summary of this function goes here
%   Detailed explanation goes here
NNE = mesh.NNE;

map = 1:3;
X = [];
Y = [];
Z = [];
C = [];
Width = max(mesh.coor(:,2));
for i=1:NNE
   node = mesh.connec(i,map);
   T = mesh.coor(node',:);
   
   u = U(3*(i-1)+(1:3),:);
   X = [X T(:,1)+alp_def.*u(:,1)];
   Y = [Y T(:,2)+alp_def.*u(:,2)];
   Z = [Z T(:,3)+alp_def.*u(:,3)];
  % C = [C Utr(3*(i-1)+(1:3))];
   C = [C u(:,3)];
end
Y_symmetry =   2*Width -Y;
colormap(jet)
patch(X,Y,Z,C,'FaceColor','interp');
patch(X,Y_symmetry,Z,C,'FaceColor','interp');

% 
% bc = mesh.bc;
% edge_u = find(bc.u==1);
% for i = 1:size(edge_u,1)
%    edg = edge_u(i);
%    pos_nd_edg = mesh.edges(:,3)==edg;
%    nd_edg = mesh.edges(pos_nd_edg,1:2);
%    xx = [mesh.coor(nd_edg(:,1),1) mesh.coor(nd_edg(:,2),1)]';
%    yy = [mesh.coor(nd_edg(:,1),2) mesh.coor(nd_edg(:,2),2)]';
%    patch(xx,yy,ones(2,size(nd_edg,1)),'FaceColor','none','EdgeColor','black','LineWidth',3.);
% end


axis equal;
X = xlim;
Y = ylim;
Z = zlim;
DX = X(2)-X(1);
DY = Y(2)-Y(1);
DZ = Z(2)-Z(1);
xlim([X(1)-0.05*DX X(2)+0.05*DX]);
ylim([Y(1)-0.05*DY Y(2)+0.05*DY]);
zlim([Z(1)-0.05*DZ Z(2)+0.05*DZ]);

end

