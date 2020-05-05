function U = pseudo_mech_cont(yd,mesh)
%PSEUDO_MECH Summary of this function goes here
%   Detailed explanation goes here
U = sparse(3*mesh.NNE,1);

NED_utr = mesh.act_edges.utrans.NED;
NED_utg = mesh.act_edges.utang.NED;
NED_un = mesh.act_edges.unorm.NED;
for e=1:mesh.NNE
    
    node = mesh.connec(e,:);
    T = mesh.coor(node(1:3),:);
    
    i1 = [1 2 3];
    i2 = [2 3 1];
    for j=1:3
        tang = T(i2(j),:)'-T(i1(j),:)';
        long = norm(tang);
        tang = tang./long;
        
       
        if (node(i1(j))>node(i2(j)));
            tang = -tang;
            temp = i1(j);
            i1(j) = i2(j);
            i2(j) = temp;
        end
        
        id_utr = mesh.act_edges.utrans.node2edg(node(i1(j)),node(i2(j)));
        id_utg = mesh.act_edges.utang.node2edg(node(i1(j)),node(i2(j)));
        id_un = mesh.act_edges.unorm.node2edg(node(i1(j)),node(i2(j)));
        
        if id_un>0
            U(3*(e-1)+[i1(j);i2(j)],1) = yd(2*(id_un-1)+(1:2))./long.*2;
        end
        if id_utg>0
            U(3*(e-1)+[i1(j);i2(j)],2) = yd(2*NED_un + 2*(id_utg-1)+(1:2))./long.*2;
        end
        if id_utr>0
            U(3*(e-1)+[i1(j);i2(j)],3) = yd(2*(NED_un+NED_utg)+2*(id_utr-1)+(1:2))./long.*2;
        end
    end

end

end

