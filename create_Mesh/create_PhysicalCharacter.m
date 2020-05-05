function textPhysic = create_PhysicalCharacter(plate_Input)

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
%% Specification of physical character

        segmentalEdge1 = num2str(2);
        for i = 3:(nJoints+2)
                segmentalEdge1 = strcat(segmentalEdge1,',',num2str(i));
        end
        
       segmentalEdge2 = num2str(nJoints+3);
        for i = (nJoints+4):(2*nJoints+3)
                segmentalEdge2 = strcat(segmentalEdge2,',',num2str(i));
        end
        
        lineListJoints = num2str(2*nJoints+4);
        for i =  (2*nJoints+5):( 3*nJoints+3)
                lineListJoints = strcat(lineListJoints,',',num2str(i));
        end
        
        surfaceList = num2str(1);
        for i = 2:nSurfaces
                surfaceList = strcat(surfaceList,',',num2str(i));
        end
        switch plate_Input.Joints.Direction
            case 'Off'
                                                                                                                                                   % Physical number
                textPhysic{1} = 'Physical Line(1) = {1};\n';                        % bottom edge of pannel                      1
                textPhysic{2} = strcat('Physical Line(2) = {2};\n');             % right edge of pannel                           2  
                textPhysic{3} = strcat('Physical Line(3) = {3};\n');              % top edge of pannel                            3
                textPhysic{4} = strcat('Physical Line(4) = {4};\n');              % left edge of pannel                             4    
                textPhysic{5} = strcat('Physical Surface(1) = {',surfaceList,'};\n'); % surfaces
            case 'Horizontal'      
                                                                                                                                                                                % Physical number
                textPhysic{1} = 'Physical Line(1) = {1};\n';                                                   % bottom edge of pannel                      1
                textPhysic{2} = strcat('Physical Line(2) = {',segmentalEdge1,'};\n');             % right edge of pannel                           2  
                textPhysic{3} = strcat('Physical Line(3) = {',num2str( 4 + 3*nJoints ),'};\n');  % top edge of pannel                            3
                textPhysic{4} = strcat('Physical Line(4) = {',segmentalEdge2,'};\n');              % left edge of pannel                             4    
                textPhysic{5} = strcat('Physical Line(5) = {',lineListJoints,'};\n');                    % joints of pannel                                5
                textPhysic{6} = strcat('Physical Surface(1) = {',surfaceList,'};\n'); % surfaces                
            case 'Vertical' 
                if nJoints ==0
                                                                                                                                                                                  % Physical number
                    textPhysic{1} = 'Physical Line(1) = {2};\n';                                                   % bottom edge of pannel                      1
                    textPhysic{2} = strcat('Physical Line(2) = {4};\n');             % right edge of pannel                           2  
                    textPhysic{3} = strcat('Physical Line(3) = {3};\n');  % top edge of pannel                            3
                    textPhysic{4} = strcat('Physical Line(4) = {1};\n');              % left edge of pannel                             4    
                    textPhysic{5} = strcat('Physical Surface(1) = {',surfaceList,'};\n'); % surfaces
                else
                    textPhysic{1} = strcat('Physical Line(1) = {',segmentalEdge1,'};\n');             % bottom edge of pannel                      1
                    textPhysic{2} = strcat('Physical Line(2) = {',num2str( 4 + 3*nJoints ),'};\n');   % right edge of pannel                         2 
                    textPhysic{3} = strcat('Physical Line(3) = {',segmentalEdge2,'};\n');               % top edge of pannel                           3 
                    textPhysic{4} =          'Physical Line(4) = {1};\n';                                          % left edge of pannel                            4  
                    textPhysic{5} = strcat('Physical Line(5) = {',lineListJoints,'};\n');                    % joints of pannel                                5
                    textPhysic{6} = strcat('Physical Surface(1) = {',surfaceList,'};\n'); % surfaces
                end
        end
            
end