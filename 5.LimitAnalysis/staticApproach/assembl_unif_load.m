function F = assembl_unif_load(mesh,pressure,orientation)
%ASSEMBL_UNIF_LOAD Assembly procedure for a load uniform/element
%   pressure = [pX;pY;pZ]; if orientation = 'global'
%   pressure = [px;py;pz]; if orientation = 'local'


F = sparse(9*mesh.NNE,1);

for e=1:mesh.NNE
   node_el = mesh.connec(e,:);
   
   T = mesh.coor(node_el,:);
   Rloc = rep_loc(T);
   
   switch orientation
       case 'global'
           orient = diag([1,1,1])*Rloc';
       case 'local'
           orient = eye(3);
       otherwise
           error('Not defined')
   end
   p = orient*pressure(e,:)';
   
   F(9*(e-1)+[1 2 9]) = p;

end
end