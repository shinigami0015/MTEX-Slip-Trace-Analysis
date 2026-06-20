function idx = match_slip_label(label, ss, cs, m)
    % MATCH_SLIP_LABEL Maps a text string slip system label back to its MTEX numeric index.
    %
    % Purpose:
    %   This helper function dynamically translates string identifiers from tables 
    %   (e.g., '(110)[1-11]' or '(10-10)[11-20]') back into their respective index 
    %   positions within the theoretical slip system array. It evaluates the 3-index 
    %   Miller notation for BCC cubic structures and the 4-index Miller-Bravais 
    %   notation for HCP hexagonal structures.
    %
    % Inputs:
    %   label - The input text string character array to match (e.g., from user CSV files)
    %   ss    - Array of theoretical MTEX slipSystem objects
    %   cs    - Cell array containing crystal symmetry objects {BCC, HCP}
    %   m     - Current active phase index indicator (1 = BCC, 2 = HCP)
    %
    % Outputs:
    %   idx   - The integer array index of the matched slip system. Returns NaN if unmatched.

    % Initialize default return index value as Not-a-Number
    idx = NaN;
    
    % Loop through every theoretical slip system variant in the active phase array
    for j = 1:length(ss)
        
        if m == 1
            % -------------------------------------------------------------
            % PHASE 1: BETA (BCC) - 3-INDEX CUBIC SYMMETRY
            % -------------------------------------------------------------
            % Convert normal vector and burger vector into standard Miller indices
            hkl = Miller(ss(j).n, cs{m}, 'plane');
            uvw = Miller(ss(j).b, cs{m}, 'direction');
            
            % Isolate and round the individual directional integers
            hkl_vals = round([hkl.h, hkl.k, hkl.l]);
            uvw_vals = round([uvw.h, uvw.k, uvw.l]);
            
            % Synthesize standard cubic string format: (hkl)[uvw]
            this_label = sprintf('(%d%d%d)[%d%d%d]', hkl_vals, uvw_vals);
            
        else
            % -------------------------------------------------------------
            % PHASE 2: ALPHA (HCP) - 4-INDEX HEXAGONAL SYMMETRY
            % -------------------------------------------------------------
            % Convert plane normal and direction into Miller-Bravais hexagonal representations
            hkl = Miller(ss(j).n, cs{m}, 'plane');
            uvw = Miller(ss(j).b, cs{m}, 'direction');
            
            % Isolate and round indices including the redundant 3rd basal index (i and T)
            hkl_vals = round([hkl.h, hkl.k, hkl.i, hkl.l]);
            uvw_vals = round([uvw.U, uvw.V, uvw.T, uvw.W]);
            
            % Synthesize hexagonal string format: (hkil)[UVTW]
            this_label = sprintf('(%d%d%d%d)[%d%d%d%d]', hkl_vals, uvw_vals);
        end
        
        % -------------------------------------------------------------
        % LABEL CORRELATION CHECK
        % -------------------------------------------------------------
        % Perform a literal string comparison between the search query and current iteration label
        if strcmp(label, this_label)
            idx = j; % Capture matching index handle position
            return;  % Immediately exit function to save compute time
        end
    end
end