function textLines   =  create_Lines(plate_Input)

%% Initialisation
    switch plate_Input.Joints.Direction
        case 'Off'
            nJoints =0;
        case 'Horizontal'
            nJoints    = plate_Input.Joints.Number;
        case 'Vertical'
            nJoints    = floor(plate_Input.Joints.Number/2);
    end
            nLines     = 4 + 3*nJoints;
            textLines = cell(nLines,1);
    %% Specification of Lines    
    textLines{1} = 'Line(1) = {1, 2};\n';
    rightNode_I = [2 5:(nJoints+4) ];
    rightNode_J = [5:(nJoints+4) 3];
    leftNode_I   = [1 (5+nJoints):(4+nJoints*2)];
    leftNode_J   = [(5+nJoints):(4+nJoints*2) 4];

    for k = 1:nJoints+1
         % line elements of one segmental edge of the rectangular plate
        textLines{k+1} = strcat('Line(',num2str(k+1),') = {',num2str(rightNode_I(k)),',',num2str(rightNode_J(k)) ,'};\n'); 
        % line elements of the other segmental edge of the rectangular plate
        textLines{k+nJoints+2} = strcat('Line(',num2str(k+nJoints+2),') = {',num2str(leftNode_I(k)),',',num2str(leftNode_J(k)) ,'};\n');
        % line elements of continuous edges 
        textLines{k+nJoints*2+3} =strcat('Line(',num2str(k+nJoints*2+3),') = {',num2str(rightNode_J(k)),',',num2str(leftNode_J(k)) ,'};\n');
    end
end