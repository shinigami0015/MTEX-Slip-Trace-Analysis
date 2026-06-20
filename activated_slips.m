function activatedInfo = activated_slips(grains_phase, grainIDs1, schmidFactors1, slipSystems1, slipFamilies1, traceAngles1)
    % ACTIVATED_SLIPS Identifies the most likely active slip system for each grain.
    % 
    % Purpose:
    %   This function loops through every grain in a given phase and identifies 
    %   the slip system that has the highest Schmid factor (the highest resolved 
    %   shear stress), which is theoretically the one most likely to activate 
    %   during mechanical testing.
    %
    % Inputs:
    %   grains_phase   - MTEX grain object containing data for the current phase
    %   grainIDs1      - List of grain IDs associated with calculated slip properties
    %   schmidFactors1 - Array of calculated Schmid factors
    %   slipSystems1   - Cell array of slip system names/labels
    %   slipFamilies1  - Cell array of slip system family classifications
    %   traceAngles1   - Array of calculated surface slip trace angles
    %
    % Outputs:
    %   activatedInfo  - A MATLAB Table summarizing the dominant slip system 
    %                    information for each grain ID.

    % Extract the unique identification numbers for all grains in this phase
    grainIDs = grains_phase.id;
    
    % Prepare output structure to temporarily hold the data for each grain
    activatedInfo = struct( ...
        'grainID', [], ...
        'schmidFactor', [], ...
        'slipSystem', [], ...
        'slipFamily', [], ...
        'traceAngle', [] ...
    );
    
    % Initialize a counter to keep track of rows in the output structure
    count = 1;
    
    % Loop through every individual grain identified in this phase
    for i = 1:length(grainIDs)
        gID = grainIDs(i);
        
        % Find all indices in the master arrays that belong to this specific grain
        idx = find(grainIDs1 == gID);
        
        % Out of all the slip systems for this grain, find the one with the maximum Schmid factor
        [~, maxIdxRel] = max(schmidFactors1(idx));
        
        % Map the relative index back to the absolute index of the master arrays
        maxIdx = idx(maxIdxRel);

        % Extract and save the parameters of the highest-ranked slip system
        activatedInfo(count).grainID      = gID;
        activatedInfo(count).schmidFactor = schmidFactors1(maxIdx);
        activatedInfo(count).slipSystem   = slipSystems1{maxIdx};
        activatedInfo(count).slipFamily   = slipFamilies1{maxIdx};
        activatedInfo(count).traceAngle   = traceAngles1(maxIdx);
        
        % Move to the next slot in our data structure
        count = count + 1;
    end
    
    % Convert the structure array into a standard MATLAB table for clean viewing and exporting
    activatedInfo = struct2table(activatedInfo);
end