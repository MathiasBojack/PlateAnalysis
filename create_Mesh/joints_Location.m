function plate_Input = joints_Location(plate_Input)
        switch plate_Input.Joints.Direction
            case 'Off'
                 plate_Input.Joints.Location =[];
            case 'Horizontal'
               % plate_Input.Joints.Number >=1
                for i = 1:plate_Input.Joints.Number
                    % vertical position of each horizontal joints
                    plate_Input.Joints.Location.Z(i) = plate_Input.Pannels.Height*i; 
                end    
            case 'Vertical'
                % plate_Input.Joints.Number >=1
                for i = 1:plate_Input.Joints.Number
                    % horizontal position of each horizontal joints
                    plate_Input.Joints.Location.Y(i) = plate_Input.Pannels.Width*i; 
                end    
        end
end