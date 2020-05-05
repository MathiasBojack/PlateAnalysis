function res = list_disc_dof(disc,mesh)
%LIST_DISC_DOF Summary of this function goes here
%   Detailed explanation goes here
nb_disc = zeros(disc.NED,1);
elem_at_disc = cell(disc.NED,1);
for id=1:disc.NED
    nodes = disc.edg2node(id,:);
    elem_at_nodes = find_el_from_node(nodes,mesh.connec);
    n = length(elem_at_nodes{1});
    nb_disc(id) = max(n-1,1);
    elem_at_disc{id} = sort(elem_at_nodes{1});
end
res.nb_disc = nb_disc;
res.elem_at_disc = elem_at_disc;
end

