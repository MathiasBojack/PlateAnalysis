function [textLineLoops,textSurfaces] = create_Surfaces(plate_Input)

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
        textLineLoops = cell(nSurfaces,1);
        textSurfaces    = cell(nSurfaces,1);
        
%% Specification of Line loops
        Line_I  = [1 -((nJoints*2+4):(nJoints*3+3)) ];  % bottom of the rectangular region
        Line_J  =  2 : (nJoints+2) ;                            % right of the rectangular region
        Line_K =  (nJoints*2+4):(nJoints*3+4) ;       % top of the rectangular region
        Line_L =  -((nJoints+3):(2*nJoints+4)) ;       % left of the rectangular region
        
        for j = 1:nSurfaces
            textLineLoops{j} = strcat('Line Loop(',num2str(j),') = {',num2str(Line_I(j)),',',num2str(Line_J(j)),',',num2str(Line_K(j)),',',num2str(Line_L(j)),'};\n');
            textSurfaces{j}    = strcat('Plane Surface(',num2str(j),') = {',num2str(j),'};\n');
        end
end