function BC = boundary_conditions(mesh)
%BOUNDARY_CONDITIONS computes the index list I such that the corresponding
%auxiliary variable mu_I has to be set to 0 because it does not contribute
%to Prm
BC = [];
% number of discontinuity terms/edge
sd = 15;

for e=1:mesh.NNE
    node = mesh.connec(e,:);
    
    i1 = [1 2 3];
    i2 = [2 3 1];
    for j=1:3
       idun = mesh.act_edges.prmun.node2edg(node(i1(j)),node(i2(j)));
       idut = mesh.act_edges.prmut.node2edg(node(i1(j)),node(i2(j)));
       iduz = mesh.act_edges.prmuz.node2edg(node(i1(j)),node(i2(j)));
       idrt = mesh.act_edges.prmrt.node2edg(node(i1(j)),node(i2(j)));
       idrn = mesh.act_edges.prmrn.node2edg(node(i1(j)),node(i2(j)));
       id = mesh.act_edges.prm.node2edg(node(i1(j)),node(i2(j)));
            
       if (id>0)
           
            elem_at_disc = mesh.id_dof.elem_at_disc{id};
            if length(elem_at_disc)==1
               id_g1 = sum(mesh.id_dof.nb_disc(1:id-1))+1; 
               
               if (idun<0)
                  BC = [BC;sd*(id_g1-1)+[1 4 7]'];
               end
               if (idut<0)
                  BC = [BC;sd*(id_g1-1)+[2 5 8]'];
               end
               if (iduz<0)
                  BC = [BC;sd*(id_g1-1)+[3 6 9]'];
               end
               if (idrn<0)
                  BC = [BC;sd*(id_g1-1)+[10;12]];
               end
               if (idrt<0)
                  BC = [BC;sd*(id_g1-1)+[11;13]];
               end
            end
            
       end
    end
    
end

end

