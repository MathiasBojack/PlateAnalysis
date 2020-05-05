function Fl = assembl_lin_load(F,line,mesh,orientation)
%ASSEMBL_LIN_LOAD Summary of this function goes here
%   Detailed explanation goes here
NED_rt = mesh.act_edges.rottan.NED;
NED_rn = mesh.act_edges.rotnorm.NED;
NED_utr = mesh.act_edges.utrans.NED;
NED_utg = mesh.act_edges.utang.NED;
NED_un = mesh.act_edges.unorm.NED;

Fl = sparse(3*(NED_rt+NED_rn)+2*(NED_utr+NED_utg+NED_un),1);

for j=1:size(line,1)
    node = line(j,:);

    T = mesh.coor(node,:);
    tang = T(2,:)'-T(1,:)';
    tang = tang./norm(tang);
    if (node(1)>node(2));
            tang = -tang;
    end
    id = mesh.act_edges.all.node2edg(node(1),node(2));
    zloc_av = mesh.edg_norm(id,:)';
    xG = mean(T);
    xM = mean(mesh.coor(node,:));
    normal_av = cross(tang,zloc_av);
    normal_av = sign((xM-xG)*normal_av).*normal_av;  
    Rloc = [normal_av tang zloc_av];
    
    switch orientation
       case 'global'
           orient = Rloc';
       case 'local'
           orient = eye(3);
       otherwise
           error('Not defined')
    end
    
    id_utr = mesh.act_edges.utrans.node2edg(node(1),node(2));
    id_utg = mesh.act_edges.utang.node2edg(node(1),node(2));
    id_un = mesh.act_edges.unorm.node2edg(node(1),node(2));
    
    FF = orient*F(2*(j-1)+(1:2),:)';
    
    if (id_un>0)
        Fl(3*(NED_rt+NED_rn) + 2*(id_un-1)+(1:2)) = FF(1,:);
    end
    if (id_utg>0)
        Fl(3*(NED_rt+NED_rn) + 2*NED_un + 2*(id_utg-1)+(1:2)) = FF(2,:);
    end    
    if (id_utr>0)
        Fl(3*(NED_rt+NED_rn) + 2*(NED_un+NED_utg) + 2*(id_utr-1)+(1:2)) = FF(3,:);
    end    
end

