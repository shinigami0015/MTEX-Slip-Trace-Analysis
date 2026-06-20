function [u, v, angle_deg, angle_deg_fin] = angle_with_horizontal(slipPlane)
    % ANGLE_WITH_HORIZONTAL Calculates the surface trace line of a slip plane.
    %
    % Purpose:
    %   This function determines where a 3D crystalline slip plane intersects 
    %   the 2D observation surface (the sample surface). It calculates the resulting 
    %   surface trace vector (u, v), finds its mathematical angle, and applies 
    %   coordinate system transformations to align it with physical microscope images.
    %
    % Inputs:
    %   slipPlane     - An MTEX vector3d object representing the slip plane normal
    %
    % Outputs:
    %   u, v          - Scaled 2D component values representing the trace line direction
    %   angle_deg     - The raw mathematical angle of the trace line (0 to 360 degrees)
    %   angle_deg_fin - The corrected, physical trace angle mapped between -90 and 90 degrees

    % Define the out-of-plane direction (Z-axis) perpendicular to the sample surface
    z = vector3d.Z;  
    
    % The intersection line (slip trace) is perpendicular to both the plane normal and the Z-axis
    slipTrace = cross(slipPlane, z);  
    
    % Project the 3D trace vector into a 2D coordinate system based on the camera view
    if z == vector3d.X
        u = slipTrace.y; v = slipTrace.z;
    elseif z == vector3d.Y
        u = slipTrace.x; v = slipTrace.z;
    else
        u = slipTrace.x; v = slipTrace.y;
    end

    % Calculate the length (magnitude) of the 2D projected trace vector
    mag = norm([u v]);
    
    % If the vector has a valid length, normalize it and scale it to a fixed length of 2.5 for plotting
    if mag > 1e-3
        u = u / mag * 2.5;
        v = v / mag * 2.5;
    else 
        % Handle edge cases where the plane is perfectly parallel to the surface
        u = 0;
        v = 0;
    end
    
    % Calculate the raw mathematical angle of the vector in degrees from the horizontal axis
    angle_deg = mod(atan2d(v, u), 360);  

    % Invert the angle direction to account for camera rotation and image flipping adjustments
    angle_deg_rot_flip = -1 * angle_deg; 
    
    % Standardize the angle so it falls strictly within the standard physical range of -90 to +90 degrees
    if angle_deg_rot_flip > 90
        angle_deg_fin = angle_deg_rot_flip - 180;
    elseif angle_deg_rot_flip >-90
        angle_deg_fin = angle_deg_rot_flip;
    elseif angle_deg_rot_flip >-180
        angle_deg_fin = angle_deg_rot_flip + 180;
    elseif angle_deg_rot_flip >-270
        angle_deg_fin = angle_deg_rot_flip + 180;
    elseif angle_deg_rot_flip >-360
        angle_deg_fin = angle_deg_rot_flip + 360;    
    else 
        angle_deg_fin = angle_deg_rot_flip;
    end
end