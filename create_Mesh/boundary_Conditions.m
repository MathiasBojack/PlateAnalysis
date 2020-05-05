function  created_Mesh = boundary_Conditions(plate_Input,created_Mesh)
% This function specifies the boundary conditions for plates with joints
% simply supported around four sides.

        % displacement normal to the boundary
        created_Mesh.bc.unorm = [1 1 0 0];
        % displacement tangent to the boundary
        created_Mesh.bc.utang = [0 0 0 0];
        % out-of-plan displacement
        created_Mesh.bc.utrans = [1 0 1 1];
        created_Mesh.bc.rottan = created_Mesh.bc.utrans;
        switch plate_Input.Joints.Direction
            case 'Off'
                    created_Mesh.bc.rotnorm = [0 1 0 0];
            case 'Horizontal'
                    created_Mesh.bc.rotnorm = [0 1 0 0];
            case 'Vertical'
                    switch mod(plate_Input.Joints.Number,2)
                        case 1
                            created_Mesh.bc.rotnorm = [0 0 0 0];
                        case 0
                            created_Mesh.bc.rotnorm = [0 1 0 0];
                    end
        end                        
end