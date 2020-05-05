function    [textLineDiscretization,textSurfaceDiscretization ]= create_Discretization(plate_Input,mesh_Info)
%% Initialisation
        switch plate_Input.Joints.Direction
            case 'Off'
                nJoints =0;
            case 'Horizontal'
                nJoints    = plate_Input.Joints.Number;
            case 'Vertical'
                nJoints    = floor(plate_Input.Joints.Number/2);
        end
        nSurfaces = 1 + nJoints;
%% Discretization


        
        switch plate_Input.Joints.Direction
            case 'Off'
                textLineDiscretization{1} = strcat('Transfinite Line{2,3} =', num2str(mesh_Info.Division.Z),';\n'); 
                textLineDiscretization{2} = strcat('Transfinite Line{1,4} =', num2str(ceil(mesh_Info.Division.Y/2)),';\n'); 
            case 'Horizontal'
                noList1 =  num2str(2);
                for i  = [3:(nJoints+2) (nJoints+3):(2*nJoints+3)]    % element number of line elements on the left and right edge
                    noList1 =  strcat( noList1, ',',num2str(i) ); 
                end
                noList2 =  num2str(1);
                for i  = [(2*nJoints+4):(3*nJoints+4) ]  % element number of horizontal line elements 
                    noList2 =  strcat( noList2, ',',num2str(i) );
                end
                % left and right edge discretization 
                textLineDiscretization{1} = strcat('Transfinite Line{', noList1,'} =', num2str(mesh_Info.Division.Z),';\n'); 
                 % horizontal element discretization 
                textLineDiscretization{2} = strcat('Transfinite Line{', noList2,'} =', num2str(ceil(mesh_Info.Division.Y/2)),';\n');
            case 'Vertical'
                if nJoints ==0  %  plate_Input.Joints.Number ==1,
                        textLineDiscretization{1} = strcat('Transfinite Line{1,4} =', num2str(mesh_Info.Division.Z),';\n'); 
                        textLineDiscretization{2} = strcat('Transfinite Line{2,3} =', num2str(mesh_Info.Division.Y),';\n'); 
                else %  plate_Input.Joints.Number >=2,
                        noList1 =  num2str(1);
                        for i  = [(2*nJoints+4):(3*nJoints+4) ]  % element number of vertical line elements 
                            noList1 =  strcat( noList1, ',',num2str(i) );
                        end
                         textLineDiscretization{1} = strcat('Transfinite Line{', noList1,'} =', num2str(mesh_Info.Division.Z),';\n'); 
                         
                        switch mod(plate_Input.Joints.Number,2)
                            case 0 %  all the joints are equally located on one side of the plane of symmetry
                                    noList2 =  num2str(2);
                                    for i  = [3:(nJoints+1) (nJoints+3):(2*nJoints+2)]    % element number of line elements on the upper and lower edge
                                        noList2 =  strcat( noList2, ',',num2str(i) ); 
                                    end
                                    textLineDiscretization{2} = ...
                                    strcat(     'Transfinite Line{', noList2,'} =', num2str(mesh_Info.Division.Y),';\n',...
                                               'Transfinite Line{', num2str(nJoints+2),',',num2str(2*nJoints+3),'} =', num2str(ceil( mesh_Info.Division.Y/2)),';\n');
                            case 1  % one joint is located on the plane of symmetry
                                    noList2 =  num2str(2);
                                    for i  = [3:(nJoints+2) (nJoints+3):(2*nJoints+3)]    % element number of line elements on the upper and lower edge
                                        noList2 =  strcat( noList2, ',',num2str(i) ); 
                                    end                                
                                    textLineDiscretization{2} = strcat('Transfinite Line{', noList2,'} =', num2str(mesh_Info.Division.Y),';\n'); 
                        end
                        % vertical element discretization 
                       
                end
            end
        
        switch mesh_Info.Direction
            case 'Right'
                Direction1 = 'Left';
                Direction2 = 'Right';
            case 'Left'
                Direction1 = 'Right';
                Direction2 = 'Left';
            case 'AlternateRight'
                Direction1 = 'AlternateRight';
                Direction2 = 'AlternateRight';
            case 'AlternateLeft'
                Direction1 = 'AlternateLeft';
                Direction2 = 'AlternateLeft';
            case 'Alternate'
                Direction1 = 'Alternate';
                Direction2 = 'Alternate';
            otherwise 
                error('invalid choice')
        end
        textSurfaceDiscretization = cell(nSurfaces,1);
        textSurfaceDiscretization{1} = strcat('Transfinite Surface{',num2str(1),'}', Direction1,';\n'); % Surfaces
        for i  = 2:nSurfaces  % element number of surface elements 
            textSurfaceDiscretization{i} = strcat('Transfinite Surface{',num2str(i),'}', Direction2,';\n'); % Surfaces
        end
        