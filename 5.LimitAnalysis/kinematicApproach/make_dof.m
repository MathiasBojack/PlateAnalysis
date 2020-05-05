function mesh = make_dof(mesh)
%MAKE_DOF Construction of the dof tables
%   the total number of dof is mesh.NEQ = 24*NNE
% and the tables for edges which are active (act_edges.prm), contribute to Prm through 
% velocity jump (act_edges.prmu), tangential rotation jumps (act_edges.prmt), normal rotation jumps (act_edges.prmn) or
% irrespective of any BCs (act_edges.all)


edges_prmun = active_edges(mesh,1-mesh.bc.unorm);
edges_prmut = active_edges(mesh,1-mesh.bc.utang);
edges_prmuz = active_edges(mesh,1-mesh.bc.utrans);
edges_prmrt = active_edges(mesh,1-mesh.bc.rottan);
edges_prmrn = active_edges(mesh,1-mesh.bc.rotnorm);
edges_prm = active_edges(mesh,and(1-mesh.bc.unorm,and(1-mesh.bc.utang,and(1-mesh.bc.utrans,and(1-mesh.bc.rottan,1-mesh.bc.rotnorm)))));
edges_all = active_edges(mesh,0.*mesh.bc.unorm);

mesh.act_edges.prmun = edges_prmun;
mesh.act_edges.prmut = edges_prmut;
mesh.act_edges.prmuz = edges_prmuz;
mesh.act_edges.prmrt = edges_prmrt;
mesh.act_edges.prmrn = edges_prmrn;
mesh.act_edges.prm = edges_prm;
mesh.act_edges.all = edges_all;

mesh.id_dof = list_disc_dof(mesh.act_edges.prm,mesh);

end

