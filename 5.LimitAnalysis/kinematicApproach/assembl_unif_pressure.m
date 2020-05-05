function F = assembl_unif_pressure(mesh,pressure,orientation)
%ASSEMBL_UNIF_LOAD Assembly procedure for a load uniform/element
%   pressure = [pX;pY;pZ]; if orientation = 'global'
%   pressure = [px;py;pz]; if orientation = 'local'


F = sparse(mesh.Nu,1);

for e=1:mesh.NNE
   node_el = mesh.connec(e,:);
   
   T = mesh.coor(node_el,:);
   Rloc = rep_loc(T);
   
   switch orientation
       case 'global'
           orient = eye(3);
       case 'local'
           orient = Rloc;
       otherwise
           error('Not defined')
   end
   p = orient*pressure(e,:)';
   area = el_area(e,mesh);
   
   F(27*(e-1)+(10:18)) = kron([1;1;1],p./3.*area);

end
end