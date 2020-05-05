function B = strain_matrix(mesh)
%STRAIN_MATRIX : Assembling of the global strain matrix B such that
 %                   B*U = d

%   if d_e collects the 18 strain variables of element nb e, the global dof vector U is of length = Ndof*NNE
%   with d_e = {chi_xx,chi_yy,chi_xy,eps_xx^1,eps_yy^1,..,2eps_xy^3,gamma_x^1,gamma_y^1,...,gamma_y^3}

NNE = mesh.NNE;
s = mesh.ndof;
Bc = cell(NNE,1);
for e=1:NNE
    node = mesh.connec(e,:);
    T = mesh.coor(node,:);
    Bc{e} = loc_strain(T,mesh);
end
B = blkdiag(Bc{:});
end

