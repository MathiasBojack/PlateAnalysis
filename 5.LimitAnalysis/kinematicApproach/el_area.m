function area = el_area(elem,mesh)
%EL_AREA Computes the area of triangular elements
n = length(elem);
area = zeros(n,1);

for i=1:n
   el = elem(i);
   node_el = mesh.connec(el,:);
   T = mesh.coor(node_el,:);
   detJ = norm(cross(T(3,:)'-T(1,:)',T(2,:)'-T(1,:)'));          % jacobian
   area(i) = detJ/2; 
end

end

