function [Utr,Utang] = pseudo_mech(y,mesh)
%PSEUDO_MECH Summary of this function goes here
%   Detailed explanation goes here
Utr = sparse(3*mesh.NNE,1);
Utang = sparse(mesh.NNE,2);

for e=1:mesh.NNE
    ye = y(9*(e-1)+(1:9));
    
    node = mesh.connec(e,:);
    T = mesh.coor(node(1:3),:);
    Ae = norm(cross(T(2,:)-T(1,:),T(3,:)-T(1,:)))./2;
    
    um = ye(9)/Ae;
    Gu1 = mean(ye(3:2:8))/Ae*6;
    Gu2 = mean(ye(4:2:8))/Ae*6;
    
    vm = ye(1:2)./Ae;
    
    
    Rloc = rep_loc(T);
    Tloc = T*Rloc(:,1:2);
    Tg = mean(Tloc);
    
    Utr(3*(e-1)+(1:3)) = (Tloc(1:3,1)-Tg(1)).*Gu1 + (Tloc(1:3,2)-Tg(2)).*Gu2+um;
    Utang(e,:) = vm;
end

end

