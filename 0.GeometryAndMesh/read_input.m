function mesh0 = read_input(fname)

fnamegeo = strcat('.\0.GeometryAndMesh\GeomFile\',fname,'.geo');
fnamemsh = strcat('.\0.GeometryAndMesh\MeshFile\',fname,'.msh');
% fid = fopen(fnamegeo, 'r');
% cmd_order = '-order 1';
% eval(horzcat('!gmsh ',cmd_order,' ',fnamegeo,' -2'));
mmesh = load_gmsh(fnamemsh);

if ((mmesh.nbTriangles == 0)&&(mmesh.nbTriangles6 == 0))
    mesh0.el_type = 'quad';
    mesh0.NNO = mmesh.nbNod;
    mesh0.NNE = mmesh.nbQuads;
    mesh0.coor = mmesh.POS(:,1:3);
    mesh0.connec = mmesh.QUADS(1:mmesh.nbQuads,1:4);
    mesh0.el_tag = mmesh.QUADS(1:mmesh.nbQuads,5);
    pos_edg = find(mmesh.LINES(:,3)>0);
    mesh0.edges = mmesh.LINES(pos_edg,:);
elseif ((mmesh.nbQuads == 0)&&(mmesh.nbTriangles6 == 0))
    mesh0.el_type = 'tri';
    mesh0.NNO = mmesh.nbNod;
    mesh0.NNE = mmesh.nbTriangles;
    mesh0.coor = mmesh.POS(:,1:3);
    mesh0.connec = mmesh.TRIANGLES(1:mmesh.nbTriangles,1:3);
    mesh0.el_tag = mmesh.TRIANGLES(1:mmesh.nbTriangles,4);
    pos_edg = find(mmesh.LINES(:,3)>0);
    mesh0.edges = mmesh.LINES(pos_edg,:);
elseif ((mmesh.nbQuads == 0)&&(mmesh.nbTriangles == 0))
    mesh0.el_type = 'tri6';
    mesh0.NNO = mmesh.nbNod;
    mesh0.NNE = mmesh.nbTriangles6;
    mesh0.coor = mmesh.POS(:,1:3);
    mesh0.connec = mmesh.TRIANGLES6(1:mmesh.nbTriangles6,1:6);
    mesh0.el_tag = mmesh.TRIANGLES6(1:mmesh.nbTriangles6,7);
    pos_edg = find(mmesh.LINES(:,3)>0);
    mesh0.edges = mmesh.LINES(pos_edg,:);
    pos_edg3 = find(mmesh.LINES3(:,3)>0);
    mesh0.edges3 = mmesh.LINES3(pos_edg3,:);
else
   error ('Triangles and Quadrangles in Mesh') 
end



clear mmesh;

end
