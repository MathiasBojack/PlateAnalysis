function file_Writing(text,plate_Input,meshFileName,file_Control)  
%% Initialisation
        switch plate_Input.Joints.Direction
            case 'Off'
                nJoints =0;
            case 'Horizontal'
                nJoints    = plate_Input.Joints.Number;
            case 'Vertical'
                nJoints    = floor(plate_Input.Joints.Number/2);
        end % number of joint
        nPoints    = 4 + 2*nJoints;
        nLines     = 4 + 3*nJoints;
        nSurfaces = 1 + nJoints;
        nPhysics = 5 +(nJoints>=1);
        %% Writing the text to the file
        file_Path =fullfile( file_Control.root_Path,'geometry_and_mesh');
        cd(file_Path);
        fid = fopen(meshFileName, 'w+');    
        for i =1:nPoints
             fprintf(fid, text.Points{i});
        end
        fprintf(fid,'\n');
        
        for k =1:nLines
            fprintf(fid,text.Lines{k});
        end
        fprintf(fid,'\n');
        
        for j = 1:nSurfaces
            fprintf(fid,text.LineLoops{j});
            fprintf(fid,text.Surfaces{j});
        end
        fprintf(fid,'\n');
        
        for l =1:nPhysics
            fprintf(fid,text.Physic{l});
        end
        fprintf(fid,'\n');
        
        for p =1:2
            fprintf(fid,text.Discretization.Line{p});
        end
        fprintf(fid,'\n');
        
        for i =1:nSurfaces
            fprintf(fid,text.Discretization.Surface{i});
        end
        fclose(fid);