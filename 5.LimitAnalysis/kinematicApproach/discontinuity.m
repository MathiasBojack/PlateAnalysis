function [Disc,edg_lgth,normal_loc] = discontinuity(mesh)
% DISCONTINUITY Matrix Disc expressing jumps of u,v,w, beta_n and beta_t
% through edges such that at a given edge : (u = normal component, v =
% tangential component, w = transversal component)
%   {u^1,v^1,...,v^3,w^3,beta_n^1,beta_t^1,beta_n^2,beta_t^2} = Disc*U

NED = mesh.NED;

s = 27;
sd = 15;

Disc = sparse(sd*NED,s*mesh.NNE);
edg_lgth = sparse(NED,1);
sign_id = ones(NED,1);
orient_id = zeros(NED,1);
normal_id = zeros(NED,3);
normal_loc = zeros(NED,2);
for e=1:mesh.NNE
    node = mesh.connec(e,:);
    Tglob = mesh.coor(node,:);
    
    Rloc = rep_loc(Tglob);
    T = Tglob*Rloc(:,1:2);
    
    i1 = [1 2 3];
    i2 = [2 3 1];
    i3 = [4 5 6];
    
    xi_eta = [1 0;0 1;0 0;1/2 1/2;0 1/2;1/2 0];
    
    for j=1:3
        
        id = mesh.act_edges.prm.node2edg(node(i1(j)),node(i2(j)));
        if id>0

            tang = Tglob(i2(j),:)'-Tglob(i1(j),:)';
            long = norm(tang);
            tang = tang./long;

            % Convention for node nb 1 and node nb 2
            if (node(i1(j))>node(i2(j)));
                tang = -tang;
                temp = i1(j);
                i1(j) = i2(j);
                i2(j) = temp;
            end

            zloc = mesh.edg_norm(id,:)';
            normal = cross(tang,zloc);
            tang = cross(zloc,normal);
            
%            rotation global to edge local frame
            R1 = [normal';tang';zloc'];
            R2 = R1;
            R3 = R1;
% 
%             zloc1 = mesh.node_norm(node(i1(j)),:)';
%             zloc2 = mesh.node_norm(node(i2(j)),:)';
%              zloc3 = (zloc1+zloc2)./norm(zloc1+zloc2);
%             normal1 = cross(tang,zloc1);
%             normal2 = cross(tang,zloc2);
%             normal3 = cross(tang,zloc3);
%             normal = mean([normal1';normal2';normal3'])';
%             tang1 = cross(zloc1,normal1);
%             tang2 = cross(zloc2,normal2);
%             tang3 = cross(zloc3,normal3);
% %            rotation global to edge local frame
%             R1 = [normal1';tang1';zloc1'];
%             R2 = [normal2';tang2';zloc2'];
%             R3 = [normal3';tang3';zloc3'];
            
            
            Proj = sparse(sd,s);
            Proj(1:3,3*(i1(j)-1)+(1:3)) = R1;
            Proj(4:6,3*(i2(j)-1)+(1:3)) = R2;
            Proj(7:9,3*(i3(j)-1)+(1:3)) = R3;
            
            switch mesh.stabilization
                case 'drilling-on'
                    xi = xi_eta(i3(j),1);
                    eta = xi_eta(i3(j),2);
                    Nthetaz = [(sqrt(2)-1)*(xi^3-xi^2) -(eta^3-eta^2) ((1-xi-eta)^3-(1-xi-eta)^2);...
                                    (xi^3-xi^2) -(sqrt(2)-1)*(eta^3-eta^2) -(1-xi-eta)^3+(1-xi-eta)^2];
                    Nthetaz2 = repmat([xi*(2*eta+xi-1) eta*(1-xi-2*eta) (1-xi-eta)*(xi-eta)],2,1);            
                    Proj(7:8,18+(3:3:9)) = R(1:2,1:2)*Nthetaz2;
                    
            end
            Proj(10:11,18+3*(i1(j)-1)+(1:3)) = [0 -1;1 0]*R1(1:2,:);
            Proj(12:13,18+3*(i2(j)-1)+(1:3)) = [0 -1;1 0]*R2(1:2,:);   
            
            Proj(14,18+3*(i1(j)-1)+(1:3)) = R1(3,:);     
            Proj(15,18+3*(i2(j)-1)+(1:3)) = R2(3,:);   
  
            
            
            Tm = mean(Tglob([i1(j);i2(j)],:));
            Tg = mean(Tglob);
            orient = sign((Tg-Tm)*normal);
            
            
            
            le = s*(e-1)+(1:s);
            
            elem_at_disc = mesh.id_dof.elem_at_disc{id};
            if length(elem_at_disc)==1
                pe1 = 1;
            else
                pe1 = find(elem_at_disc(1:end-1)==e);
            end
            
            if not(isempty(pe1))
                id_g1 = sum(mesh.id_dof.nb_disc(1:id-1))+pe1;
                Disc(sd*(id_g1-1)+(1:sd),le) = Disc(sd*(id_g1-1)+(1:sd),le) + orient.*Proj; 
                sign_id(id_g1)=-1;
                edg_lgth(id_g1) = long;
                
                if (orient_id(id_g1)==0)
                    orient_id(id_g1) = orient;
                end
            normal_id(id_g1,:) = normal_id(id_g1,:)+orient;
                
            end

            if length(elem_at_disc)>1
                pe2 = find(elem_at_disc(2:end)==e);
                if not(isempty(pe2))
                    id_g2 = sum(mesh.id_dof.nb_disc(1:id-1))+pe2;
                    Disc(sd*(id_g2-1)+(1:sd),le) = Disc(sd*(id_g2-1)+(1:sd),le) + orient.*Proj; 
                    sign_id(id_g2)=-1;
                    edg_lgth(id_g2) = long;
                
                    if (orient_id(id_g2)==0)
                        orient_id(id_g2) = orient;
                    end
                    normal_id(id_g2,:) = normal_id(id_g2,:)+orient;
                end
            end
            
            normal_loc(id,:) = normal'*Rloc(:,1:2);
        end
    end
    
end
% Disc = kron(speye(15),sparse(diag(orient_id)))*Disc;
end

