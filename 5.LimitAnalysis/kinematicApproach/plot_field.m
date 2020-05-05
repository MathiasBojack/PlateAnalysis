function res = plot_field(V,U,alp_def,mesh)
%PLOT_MESH Summary of this function goes here
%   Detailed explanation goes here

hold on
X = [];
Y = [];
Z = [];
C = [];
for e=1:mesh.NNE
   node = mesh.connec(e,:);
   T = mesh.coor(node',:);
   
   if (size(U,1)==18*mesh.NNE)
       ux = U(18*(e-1)+[1 4 7]);
       uy = U(18*(e-1)+[2 5 8]);
       uz = U(18*(e-1)+[3 6 9]);
   elseif (size(U,1)==9*mesh.NNE)
       ux = U(9*(e-1)+[1 4 7]);
       uy = U(9*(e-1)+[2 5 8]);
       uz = U(9*(e-1)+[3 6 9]);
   else
       error('Wrong dimension in U');
   end
   
   if (size(V,1) == 6*mesh.NNE)
       v = V(6*(e-1)+(1:3));
   elseif (size(V,1) == 3*mesh.NNE)
       v = V(3*(e-1)+(1:3));
   elseif (size(V,1) == mesh.NNE)
       v = repmat(V(e),3,1);
   else
       error('Wrong dimension in V')
   end
   
   X = [X T(:,1) + alp_def.*ux];
   Y = [Y T(:,2) + alp_def.*uy];
   Z = [Z T(:,3) + alp_def.*uz];

   C = [C  v];
   
end
colormap(jet)
patch(X,Y,Z,C,'FaceColor','interp','CData',C,'EdgeColor','black','LineWidth',0.2,'CDataMapping','scaled')
axis equal
end

