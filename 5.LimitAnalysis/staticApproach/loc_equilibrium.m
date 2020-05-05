function He = loc_equilibrium(Tglob)
%LOC_EQUILIBRIUM local equations of equilibrium for the quadratic
%moment/linear shear static finite element
% div M + V = 0 at the three nodes of the triangle (6 equations)
% div V = p     in the element (1 equation)
%   Detailed explanation goes here
Rloc = rep_loc(Tglob);
T =Tglob*Rloc(:,1:2);

J = [(T(1,:)-T(3,:))' (T(2,:)-T(3,:))'];   % jacobian matrix
detJ=J(1,1)*J(2,2)-J(1,2)*J(2,1);          % jacobian
invJ=1/detJ*[ J(2,2) -J(1,2); ...          % inverse jacobian matrix
       -J(2,1)  J(1,1)];
   
DNv = [1 0 -1;0 1 -1];
GN = invJ'*DNv;
div_V = [GN(1,1) GN(2,1) GN(1,2) GN(2,2) GN(1,3) GN(2,3)];

div_N = [GN(1,1) 0 GN(2,1) GN(1,2) 0 GN(2,2) GN(1,3) 0 GN(2,3);
            0 GN(2,1) GN(1,1) 0 GN(2,2) GN(1,2) 0 GN(2,3) GN(1,3)]; 

% loop on triangle end nodes
loc_coor = [1 0;0 1;0 0];
div_M_node = [];
for i=1:3
    xi = loc_coor(i,1);
    eta = loc_coor(i,2);
    
    [~,DNm] = quadratic_shape_fct(xi,eta);
    
    GN = invJ'*DNm';
    div_M = sparse(2,18);
    div_M(1,1:3:end) = GN(1,:);
    div_M(1,3:3:end) = GN(2,:);
    div_M(2,3:3:end) = GN(1,:);
    div_M(2,2:3:end) = GN(2,:);
    
    div_M_node = [div_M_node;div_M];
end

He = [div_N sparse(2,24);sparse(6,9) div_M_node speye(6);sparse(1,27) div_V];

end

