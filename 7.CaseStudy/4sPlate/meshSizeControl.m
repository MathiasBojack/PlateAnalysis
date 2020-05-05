numDiv = 20;
lc = a/numDiv;

fname = [prjRoot '\0.GeometryAndMesh\GeomFile\Geometry_a_' num2str(a) '_b_' num2str(b) '.geo' ];
fname_temp = [prjRoot '\0.GeometryAndMesh\GeomFile\Geometry_a_' num2str(a) '_b_' num2str(b) 'temp.geo' ];

fid = fopen(fname,'r+');
fid_temp = fopen(fname_temp ,'w+');

while ~feof(fid)
    tline = fgetl(fid);
    if contains(tline, 'lc=') 
        text = ['lc=' num2str(lc) ';'];
        fprintf(fid_temp, '%s', text); fprintf(fid_temp, '\n');
    else 
        fprintf(fid_temp, '%s', tline); fprintf(fid_temp, '\n');
    end

end

fclose(fid);
fclose(fid_temp);
fclose('all');
movefile(fname_temp,fname)

%% generate the mesh file

cmd_txt = ['D:\Gmsh\gmsh -order 1 ',fname, ' -2'];
system(cmd_txt)

mshName = [prjRoot '\0.GeometryAndMesh\GeomFile\Geometry_a_' num2str(a) '_b_' num2str(b) '.msh' ];
mshNameMoved = [prjRoot '\0.GeometryAndMesh\MeshFile\Geometry_a_' num2str(a) '_b_' num2str(b) '.msh' ];
movefile(mshName,mshNameMoved)

