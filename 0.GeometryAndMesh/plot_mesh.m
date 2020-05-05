function [ output_args ] = plot_mesh(mesh)
%PLOT_MESH Summary of this function goes here
%   Detailed explanation goes here
NNE = mesh.NNE;
clf
hold on

map = 1:3;
X = [];
Y = [];
Z = [];
C = [];
for i=1:NNE
   node = mesh.connec(i,map);
   T = mesh.coor(node',:);
   Rloc = rep_loc(T);
   xG = mean(T,1);
   lg_el = sqrt(norm(cross(T(2,:)-T(1,:),T(3,:)-T(1,:)))/2);
   normal = Rloc(:,3).*lg_el/2;
   a1 = Rloc(:,1).*lg_el/2;
   a2 = Rloc(:,2).*lg_el/2;
   
%    plot3([xG(1);xG(1)+normal(1)],[xG(2);xG(2)+normal(2)],[xG(3);xG(3)+normal(3)],'-b','LineWidth',2)
%    plot3([xG(1);xG(1)+a1(1)],[xG(2);xG(2)+a1(2)],[xG(3);xG(3)+a1(3)],'-r','LineWidth',2)
%    plot3([xG(1);xG(1)+a2(1)],[xG(2);xG(2)+a2(2)],[xG(3);xG(3)+a2(3)],'-g','LineWidth',2)
   
   X = [X T(:,1)];
   Y = [Y T(:,2)];
   Z = [Z T(:,3)];
   vertex = size(T,1);
   C = [C ones(vertex,1)];
end
patch(X,Y,Z,C,'FaceColor','none');


nd_edg = mesh.edges(:,1:2);
xx = [mesh.coor(nd_edg(:,1),1) mesh.coor(nd_edg(:,2),1)]';
yy = [mesh.coor(nd_edg(:,1),2) mesh.coor(nd_edg(:,2),2)]';
zz = [mesh.coor(nd_edg(:,1),3) mesh.coor(nd_edg(:,2),3)]';
patch(xx,yy,zz,ones(2,size(nd_edg,1)),'FaceColor','none','EdgeColor','black','LineWidth',2.);



view(3)
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

