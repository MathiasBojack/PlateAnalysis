function textPoints = create_Points(plate_Input)

%% Initialisation
    b               = plate_Input.Geometries.b;
    a               =  plate_Input.Geometries.a;   
        %% Horizontal joints
    switch plate_Input.Joints.Direction
        case 'Off'
                textPoints{1} = 'Point(1) = {0, 0, 0};\n';
                textPoints{2} = strcat( 'Point(2) = {0,',num2str(b/2),', 0};\n');
                textPoints{3} = strcat( 'Point(3) = {0,',num2str(b/2),',',num2str(a),'};\n');
                textPoints{4} = strcat( 'Point(4) = {0,0,',num2str(a),'};\n');
        case 'Horizontal'
                nJoints       = plate_Input.Joints.Number;
                nPoints      = 4 + 2*nJoints;
                textPoints  = cell(nPoints,1);
                %% Specification of points
                textPoints{1} = 'Point(1) = {0, 0, 0};\n';
                textPoints{2} = strcat( 'Point(2) = {0,',num2str(b/2),', 0};\n');
                textPoints{3} = strcat( 'Point(3) = {0,',num2str(b/2),',',num2str(a),'};\n');
                textPoints{4} = strcat( 'Point(4) = {0,0,',num2str(a),'};\n');

                for i =1 : nJoints
                    textPoints{i+4} = strcat( 'Point(' , num2str(i+4) , ') = {0,' ,num2str(b/2), ',' , num2str(plate_Input.Joints.Location.Z(i)) ,'};\n');
                    textPoints{i+4+nJoints} = strcat('Point(',num2str(i+4+nJoints),') = {0,0,', num2str(plate_Input.Joints.Location.Z(i)) ,'};\n');
                end 
    %% Vertical Joints  
        case 'Vertical'
            nJoints       = floor(plate_Input.Joints.Number/2);
            nPoints      = 4 + 2*nJoints;
            textPoints  = cell(nPoints,1);
            %% Specification of points
            textPoints{1} = strcat( 'Point(1) = {0,0,',num2str(a),'};\n');
            textPoints{2} = 'Point(2) = {0, 0, 0};\n';
            textPoints{3} = strcat( 'Point(3) = {0,',num2str(b/2),', 0};\n');
            textPoints{4} = strcat( 'Point(4) = {0,',num2str(b/2),',',num2str(a),'};\n');
            for i =1 : nJoints
                textPoints{i+4} = strcat('Point(',num2str(i+4),') = {0,', num2str(plate_Input.Joints.Location.Y(i)) ,',0};\n');
                textPoints{i+4+nJoints} = strcat('Point(',num2str(i+4+nJoints),') = {0,', num2str(plate_Input.Joints.Location.Y(i)) ,',',num2str(a),'};\n');
            end 
    end
end