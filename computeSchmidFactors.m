function schmidFactors = computeSchmidFactors(orientations, stressDirection, slipSystems)
    % COMPUTESCHMIDFACTORS Calculates the Schmid factor for every grain and slip system.
    %
    % Purpose:
    %   This function determines the Schmid factor ($m = \cos\phi \cdot \cos\lambda$) 
    %   for all available grains. It rotates the theoretical slip directions (b) 
    %   and slip plane normals (n) from the crystal frame to the sample frame based 
    %   on each grain's orientation, then calculates their geometric relationship 
    %   relative to the applied external loading axis (stress direction).
    %
    % Inputs:
    %   orientations    - MTEX orientation or rotation object for all grains
    %   stressDirection - A vector3d object defining the external tensile/compressive axis
    %   slipSystems     - An array of MTEX slipSystem objects
    %
    % Outputs:
    %   schmidFactors   - A matrix of size (n_grains x n_slipSystems) storing the 
    %                     resolved shear stress factor (ranging from 0 to 0.5)

    % Get the total number of grains and total number of slip systems to evaluate
    n_grains = length(orientations);
    n_ss = length(slipSystems);

    % Pre-allocate the matrix with zeros to optimize speed and memory usage
    schmidFactors = zeros(n_grains, n_ss);

    % Convert orientation data into a formal MTEX rotation matrix object
    oriMat = rotation(orientations);  

    % Loop through each individual slip system
    for j = 1:n_ss

        % --- Explicitly ensure the variables are vector3d objects ---
        % Convert the slip direction (b) and slip plane normal (n) to vector3d objects 
        % to prevent syntax errors, then rotate them into the sample reference frame.
        b_rot = rotate(vector3d(slipSystems(j).b), oriMat);  % Rotated slip direction
        n_rot = rotate(vector3d(slipSystems(j).n), oriMat);  % Rotated slip plane normal

        % --- Compute Schmid Factor ---
        % Formula: m = | (b_rotated · stress_axis) * (n_rotated · stress_axis) |
        % This uses vectorized dot products across all grains simultaneously for speed.
        schmidFactors(:, j) = abs(dot(normalize(b_rot), normalize(stressDirection)) .* dot(normalize(n_rot), normalize(stressDirection)));

    end
end