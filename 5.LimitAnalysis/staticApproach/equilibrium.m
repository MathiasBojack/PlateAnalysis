function H = equilibrium(mesh)
%EQUILIBRIUM Summary of this function goes here
%   Detailed explanation goes here

NNE = mesh.NNE;
Hc = cell(NNE,1);

for e=1:NNE
    node = mesh.connec(e,:);
    T = mesh.coor(node,:);
    Hc{e} = loc_equilibrium(T);
end

H = blkdiag(Hc{:});
end

