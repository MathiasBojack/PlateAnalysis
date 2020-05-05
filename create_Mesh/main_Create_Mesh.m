function [meshFileName,plate_Input]= main_Create_Mesh( plate_Input,mesh_Info,file_Control)
% This function create the mesh of plate by using Gmsh
% Input list:
%   plate_Input 
%           Geometries
%                   Width  : a vecteur of the width of pannels 
%                   Height:  a vecteur of the heigth of pannels
%                   Thickness
%               Pannel
%                       Height
%                       Width
%                       Center
%                       Upper_Boundary
%                       Lower_Boundary
%           Joints
%                   Direction = [ 'Horizontal', 'Vertical']
%                   Number
%                   Location
%   mesh_Info
%           Division
%                   Y       = interger number  % horizontal discretization
%                   Z        = interger number  % vertical discretization
%           Direction = [ Right, Left, Alternate,AlternateRight,AlternateLeft]
%
%  file_Control.
%           .fileName: file name of the mesh  
%                   .geo  geometry
%                   .msh  mesh
%           .root_Path 

%% Filecontrol
if isfield(file_Control,'fileName')==0
    file_Control.fileName = 'Geometry';
end

meshFileName.geo = strcat(file_Control.fileName,'.geo');
meshFileName.msh = strcat(file_Control.fileName,'.msh');

%%  Pannel with joints
% if plate_Input.Joints.Number >=1
        %% Specification of Points
        text.Points  =  create_Points(plate_Input);
        %% Specification of Lines
        text.Lines   =  create_Lines(plate_Input);
        %% Specification of Surfaces
        [text.LineLoops,text.Surfaces] = create_Surfaces(plate_Input);
        %% Specification of PhysicalCharacters
        text.Physic = create_PhysicalCharacter(plate_Input);
        %% Specification of Discretization
        [text.Discretization.Line,text.Discretization.Surface ]= create_Discretization(plate_Input,mesh_Info);
        %% 
        file_Writing(text,plate_Input,meshFileName.geo,file_Control);
%%  Pannel without joints
% elseif plate_Input.Joints.Number ==0
%         %% Initialization
%         a          = plate_Input.Geometries.a;
%         b         = plate_Input.Geometries.b;
%         %% Specification of Points
%     
%         text.Points{1} = 'Point(1) = {0, 0, 0};\n';
%         text.Points{2} = strcat( 'Point(2) = {0,',num2str(a),', 0};\n');
%         text.Points{3} = strcat( 'Point(3) = {0,',num2str(a),',',num2str(b),'};\n');
%         text.Points{4} = strcat( 'Point(4) = {0,0,',num2str(b),'};\n');
% 
%         %% Specification of Lines
%         node_I = 1:4;
%         node_J = [2:4 1];
%         for k = 1:4
%         	text.Lines{k} = strcat('Line(',num2str(k),') = {',num2str(node_I(k)),',',num2str(node_J(k)) ,'};\n');  
%         end
%         
%         %% Specification of Surfaces
%         
%         text.LineLoops{1} = 'Line Loop(1) = {1,2,3,4};\n';
%         text.Surfaces{1}    = 'Plane Surface(1) = {1};\n';
%         
%         %% Specification of PhysicalCharacters
%         text.Physic{1} =  'Physical Line(1) = {1};\n';
%         text.Physic{2} =  'Physical Line(2) = {2};\n';
%         text.Physic{3} =  'Physical Line(3) = {3};\n';
%         text.Physic{4} =  'Physical Line(4) = {4};\n';
%         text.Physic{5} =  'Physical Surface(1) = {1};\n';
% 
%         %% Specification of Discretization
%         text.Discretization.Line{1} =  strcat('Transfinite Line{1,3} =', num2str(mesh_Info.Division.Y),';\n'); 
%         text.Discretization.Line{2} =  strcat('Transfinite Line{2,4} =', num2str(mesh_Info.Division.Z),';\n'); 
%         text.Discretization.Surface{1} = strcat('Transfinite Surface{1}', mesh_Info.Direction,';\n');
%         %% 
%         file_Writing(text,plate_Input,meshFileName.geo,file_Control);
% end        
%% execution
cmd_order = '-order 1';
eval(horzcat('!gmsh ',cmd_order,' ',meshFileName.geo,' -2'));
end

