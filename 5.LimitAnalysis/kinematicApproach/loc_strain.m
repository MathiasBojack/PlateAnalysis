function Be = loc_strain(Tglob,mesh)
%LOC_STRAIN local equations of strains for the quadratic
% displacement/linear rotation kinematic finite element
% chi = \nabla_s \beta in the element (3 equations)
% eps = \nabla_s ui*a_i           at the three nodes of the triangle (9 equations)
% \gamma = \nabla w - \beta   at the three nodes of the triangle (6 equations)
%   collected as :
%   [chi_xx,chi_yy,2*chi_xy,eps_xx^1,eps_yy^1,..,2eps_xy^3,gamma_x^1,gamma_y^1,...,gamma_y^3]

Rloc = rep_loc(Tglob);
Rot = rot_loc2glob(Rloc);

T =Tglob*Rloc(:,1:2);

J = [(T(1,:)-T(3,:))' (T(2,:)-T(3,:))'];   % jacobian matrix
detJ=J(1,1)*J(2,2)-J(1,2)*J(2,1);          % jacobian
invJ=1/detJ*[ J(2,2) -J(1,2); ...          % inverse jacobian matrix
       -J(2,1)  J(1,1)];
   
DNb = [1 0 -1;0 1 -1];
GN = invJ'*DNb;
curv = sparse([GN(1,1) 0 0 GN(1,2) 0 0 GN(1,3) 0 0; ...
                     0 GN(2,1) 0 0 GN(2,2) 0 0 GN(2,3) 0; ...
                     [GN(2,1) GN(1,1) 0 GN(2,2) GN(1,2) 0 GN(2,3) GN(1,3) 0]])*kron(speye(3),[0 -1 0;1 0 0;0 0 1]);

% loop on triangle end nodes
loc_coor = [1 0;0 1;0 0];
memb = [];
memb_thetaz = [];
for i=1:3
    xi = loc_coor(i,1);
    eta = loc_coor(i,2);

    [~,DNT6] = quadratic_shape_fct(xi,eta);

    GN = invJ'*DNT6';
    
    DNthetazu = [(sqrt(2)-1)*(3*xi^2-2*xi) 0 -(3*(1-xi-eta)^2-2*(1-xi-eta))/1;...
                        0 -(3*eta^2-2*eta) -(3*(1-xi-eta)^2-2*(1-xi-eta))/1];
    DNthetazv = [(3*xi^2-2*xi) 0 (3*(1-xi-eta)^2-2*(1-xi-eta))/1;...
                        0 -(sqrt(2)-1)*(3*eta^2-2*eta) (3*(1-xi-eta)^2-2*(1-xi-eta))/1];
                    
    DNthetaz2 = [2*eta+2*xi -eta 1-2*xi;...
                       2*xi -4*eta+1-xi -1+2*eta];
                    
    GNthetazu = invJ'*DNthetaz2;                
    GNthetazv = invJ'*DNthetaz2;                
                    
    memb = [memb;kron(GN(1,:),[1 0 0]);...
             kron(GN(2,:),[0 1 0]);...
             kron(GN(1,:),[0 1 0]) + kron(GN(2,:),[1 0 0])];
         
    memb_thetaz = [memb_thetaz; kron(GNthetazu(1,:),[0 0 1]);...
                           kron(GNthetazv(2,:),[0 0 1]);...
                           kron(GNthetazu(2,:)+GNthetazv(1,:),[0 0 1])];
              
end         
                       
loc_coor = [1 0;0 1;0 0];
% loc_coor = [2/3 1/6;1/6 2/3;1/6 1/6];
% loc_coor = [1/2 1/2;0 1/2;1/2 0];
nabla_w = [];
for i=1:3
    xi = loc_coor(i,1);
    eta = loc_coor(i,2);

    [~,DNT6] = quadratic_shape_fct(xi,eta);

    GN = invJ'*DNT6';    
    nabla_w = [nabla_w;kron(GN(1,:),[0 0 1]);...
                    kron(GN(2,:),[0 0 1])];
end



switch mesh.stabilization
    case 'drilling-off'
       Be = [sparse(3,18) curv;...
        memb sparse(9,9);...
        nabla_w -kron(speye(3),[0 -1 0;1 0 0])]*Rot';
        
    case 'drilling-on'
       Be = [sparse(3,18) curv;...
        memb memb_thetaz;...
        nabla_w -kron(speye(3),[0 -1 0;1 0 0])]*Rot';
end


end

