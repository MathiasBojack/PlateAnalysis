function Liais = liaison(mesh,id_list)
%CONTINUITY Summary of this function goes here
%   Detailed explanation goes here

NEQ = mesh.NEQ;
Liais = sparse(3*length(id_list),NEQ);

for e=1:mesh.NNE
    node = mesh.connec(e,:);
    T = mesh.coor(node,:);
    
    i1 = [1 2 3];
    i2 = [2 3 1];
    i3 = [4 5 6];
    for j=1:3
        tang = T(i2(j),:)'-T(i1(j),:)';
        tang = tang./norm(tang);
        
        Rloc = rep_loc(T);
        a1 = Rloc(:,1);
        a2 = Rloc(:,2);
        zloc = Rloc(:,3);
        
        if (node(i1(j))>node(i2(j)));
            tang = -tang;
            temp = i1(j);
            i1(j) = i2(j);
            i2(j) = temp;
        end
        
        id = mesh.act_edges.all.node2edg(node(i1(j)),node(i2(j)));
        zloc_av = mesh.edg_norm(id,:)';
        normal_av = cross(tang,zloc_av);
        normal = cross(tang,zloc);
        
        if ismember(id,id_list)
        
            % expression of M_{nn},M_{nt} from their components M_{ab} in the
            % local basis of the element
            Mnn = [(normal'*a1)*(normal_av'*a1) (normal'*a2)*(normal_av'*a2) (normal'*a1)*(normal_av'*a2)+(normal_av'*a1)*(normal'*a2)];
            %Mnn = [(normal'*a1)^2 (normal'*a2)^2 2*(normal'*a1)*(normal'*a2)];
            Mnt = [(normal'*a1)*(tang'*a1) (normal'*a2)*(tang'*a2) (normal'*a1)*(tang'*a2)+(normal'*a2)*(tang'*a1)];
            
            ProjM = sparse(3,18);
            ProjM(1,3*(i1(j)-1)+(1:3)) = Mnn;
            ProjM(2,3*(i2(j)-1)+(1:3)) = Mnn;
            ProjM(3,3*(i3(j)-1)+(1:3)) = Mnn;

        
            nu1 = normal'*a1;
            nu1_av = normal_av'*a1;
            t1 = tang'*a1;
            z1 = zloc'*a1;
            z1_av = zloc_av'*a1;
            nu2 = normal'*a2;
            nu2_av = normal_av'*a2;
            t2 = tang'*a2;
            z2 = zloc'*a2;
            z2_av = zloc_av'*a2;
            ProjR = sparse(6,15);
            l1 = [3*(i1(j)-1)+(1:3) 9+2*(i1(j)-1)+(1:2)];
            ProjR(1,l1) = [nu1*nu1_av nu2*nu2_av nu2*nu1_av+nu2_av*nu1 (zloc'*normal_av)*nu1 (zloc'*normal_av)*nu2];
            ProjR(2,l1) = [nu1_av*t1 nu2_av*t2 nu2_av*t1+nu1_av*t2 (zloc_av'*tang)*nu1_av (zloc_av'*tang)*nu2_av];
            ProjR(3,l1) = [nu1*z1_av nu2*z2_av nu2*z1_av+nu1*z2_av nu1*(zloc'*zloc_av) nu2*(zloc'*zloc_av)];
            l2 = [3*(i2(j)-1)+(1:3) 9+2*(i2(j)-1)+(1:2)];
            ProjR(4:6,l2) = ProjR(1:3,l1);
        
            leM = 33*(e-1)+9+(1:18);
            leNV = 33*(e-1)+ [1:9 28:33];
        
            L = sparse(3,NEQ);
            L(:,leM) = ProjM(1:3,:); 
            LL = sparse(2,NEQ);
            LL(:,leNV) = ProjR([2 5],:);
            Liais =  [Liais;L;LL];
        end
          id_list = setdiff(id_list,id);
    end
    
end

end

