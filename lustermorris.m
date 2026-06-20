function [mprime_segments, gB, mprime_table] = lustermorris(grains, activatedInfo, cs, m)
    % LUSTERMORRIS Calculates and visualizes the Luster-Morris slip compatibility parameter.
    %
    % Purpose:
    %   This function assesses geometric slip transmission across grain boundaries.
    %   For each neighboring grain pair, it checks their active slip systems, 
    %   transforms the crystal plane normals (n) and slip directions (b) into the 
    %   sample reference frame, and calculates the Luster-Morris factor:
    %   m' = abs(cos(psi) * cos(lambda)). High values (near 1) mean the slip systems
    %   are well-aligned across the boundary, facilitating continuous deformation.
    %
    % Inputs:
    %   grains        - MTEX grain object containing individual grain data
    %   activatedInfo - MATLAB Table containing the designated active slip system per grain
    %   cs            - Cell array containing crystal symmetry objects
    %   m             - Numerical index targeting the active crystal phase
    %
    % Outputs:
    %   mprime_segments - Vector containing m' values mapped onto individual boundary lines
    %   gB              - The updated MTEX boundary object containing the geometric data
    %   mprime_table    - Summarized compilation table detailing grain-to-grain boundary values

    % 1. Extract Valid Boundaries (Discard outer scan edges and NaN regions)
    gB = grains.boundary;
    valid_mask = all(gB.grainId > 0, 2) & ~any(isnan(gB.grainId), 2);
    gB_inner = gB(valid_mask);
    
    % 2. Identify Unique Grain Pairs (Avoid duplicate calculations for shared walls)
    all_pairs = gB_inner.grainId;
    [uniquePairs, ~, idx_map] = unique(sort(all_pairs, 2), 'rows');
    
    nUnique = size(uniquePairs, 1);
    
    % Initialize storage arrays for calculation outputs
    mprime_pair_values = nan(nUnique, 1); 
    slipSystem_names_1 = strings(nUnique, 1); 
    slipSystem_names_2 = strings(nUnique, 1); 
    
    % Get full list of available theoretical slip systems
    [ss, ~] = slip_systems(cs, m);
    
    % 3. Loop over Unique Pairs
    fprintf('Calculating Luster-Morris for %d unique grain pairs...\n', nUnique);
    
    for i = 1:nUnique
        gid1 = uniquePairs(i, 1);
        gid2 = uniquePairs(i, 2);
        
        % Retrieve activated slip system info for these grains
        row1 = activatedInfo(activatedInfo.grainID == gid1, :);
        row2 = activatedInfo(activatedInfo.grainID == gid2, :);
        
        if isempty(row1) || isempty(row2)
            continue;
        end
        
        % --- STORE NAMES ---
        % Capture the names before we convert them to indices
        name1 = row1.slipSystem;
        name2 = row2.slipSystem;
        
        % Handle cell array vs string vs char data type variations safely
        if iscell(name1), name1 = string(name1{1}); end
        if iscell(name2), name2 = string(name2{1}); end
        if ischar(name1), name1 = string(name1); end
        if ischar(name2), name2 = string(name2); end
        
        % Save to our storage arrays
        slipSystem_names_1(i) = name1;
        slipSystem_names_2(i) = name2;
        
        % --- MATH AND GEOMETRY ---
        % Match the character labels to index integers in the main slip system definitions
        ss_idx1 = match_slip_label(row1.slipSystem, ss, cs, m);
        ss_idx2 = match_slip_label(row2.slipSystem, ss, cs, m);
        
        if isnan(ss_idx1) || isnan(ss_idx2)
            continue;
        end
        
        sS1 = ss(ss_idx1);
        sS2 = ss(ss_idx2);
        
        % 1. Get Mean Orientations
        ori1 = grains(gid1).meanOrientation;
        ori2 = grains(gid2).meanOrientation;
        
        % 2. Rotate vectors to Sample Frame
        n1_sample = normalize(ori1 * sS1.n);
        b1_sample = normalize(ori1 * sS1.b);
        n2_sample = normalize(ori2 * sS2.n);
        b2_sample = normalize(ori2 * sS2.b);
        
        % 3. Calculate Luster-Morris Factor
        % psi is the angle between slip plane normals; lambda is the angle between slip directions
        cos_psi = dot(n1_sample, n2_sample);
        cos_lambda = dot(b1_sample, b2_sample);
        
        mprime_pair_values(i) = abs(cos_psi * cos_lambda);
    end
    
    % 4. Map Results back to Boundary Segments for Plotting
    mprime_inner_segments = mprime_pair_values(idx_map);
    
    mprime_segments = nan(length(gB), 1);
    mprime_segments(valid_mask) = mprime_inner_segments;
    
    gB.prop.mprime = mprime_segments;
    
    % 5. Create Output Table (NOW WITH NAMES)
    valid_rows = ~isnan(mprime_pair_values);
    
    mprime_table = table(uniquePairs(valid_rows,1), uniquePairs(valid_rows,2), ...
                         mprime_pair_values(valid_rows), ...
                         slipSystem_names_1(valid_rows), ...
                         slipSystem_names_2(valid_rows), ...
        'VariableNames', {'GrainID_1', 'GrainID_2', 'LusterMorrisFactor', 'SlipSystem_1', 'SlipSystem_2'});
    
    % 6. Visualization
    figure; % Ensure a new figure is created
    colormap(hot);
    plot(gB(~isnan(mprime_segments)), 'linewidth', 2, 'property', 'mprime');
    mtexColorbar('title','Luster-Morris Factor (m'')');
    clim([0 1]); 
    hold on;
    
    % Overlay grain identification index numbers at the center centroids
    for g = 1:length(grains)
       if ~isnan(grains(g).id)
           text(grains(g).centroid.x, grains(g).centroid.y, num2str(grains(g).id), ...
               'FontSize', 8, 'Color', 'blue', 'HorizontalAlignment', 'center');
       end
    end
    hold off;
end