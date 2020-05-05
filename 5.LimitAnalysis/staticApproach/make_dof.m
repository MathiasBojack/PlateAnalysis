function mesh = make_dof(mesh)
%MAKE_DOF Summary of this function goes here
%   Detailed explanation goes here

NNE = mesh.NNE;

mesh.NEQ = 33*NNE;

edges_rottan = active_edges(mesh,mesh.bc.rottan);
edges_rotnorm = active_edges(mesh,mesh.bc.rotnorm);
edges_utrans = active_edges(mesh,mesh.bc.utrans);
edges_utang = active_edges(mesh,mesh.bc.utang);
edges_unorm = active_edges(mesh,mesh.bc.unorm);
edges_all = active_edges(mesh,0.*mesh.bc.utrans);
bc_edges = edges_on_boundary(mesh);

mesh.act_edges.rottan = edges_rottan;
mesh.act_edges.rotnorm = edges_rotnorm;
mesh.act_edges.utrans = edges_utrans;
mesh.act_edges.utang = edges_utang;
mesh.act_edges.unorm = edges_unorm;
mesh.act_edges.all = edges_all;
mesh.bc_edges = bc_edges;
end

