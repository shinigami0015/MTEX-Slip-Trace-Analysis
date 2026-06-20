function [schmidFactors, grainIDs1, schmidFactors1, slipSystems1, slipFamilies1, traceAngles1, traceAngles1_fin, q, qx, qy] = slip_trace(ebsd_phase, grains_phase, loaddir, cs, m)
    % SLIP_TRACE Calculates Schmid factors and traces active crystallographic slip planes.
    %
    % Purpose:
    %   This core module iterates over every localized grain within an active phase to
    %   compute Schmid factors, rank systems from highest to lowest resolved shear
    %   stress, and map trace projection lines onto a sample coordinate window. It
    %   automatically renders the top rank system as a solid vector line and subsequent
    %   variants as dashed line styles.
    %
    % Inputs:
    %   ebsd_phase   - Isolated phase subset of the processed MTEX EBSD dataset
    %   grains_phase - Extracted structural grain handles for the target phase
    %   loaddir      - Macroscopic load vector definition (e.g., [1, 0, 0])
    %   cs           - Cell array containing crystal symmetry objects {BCC, HCP}
    %   m            - Current active phase index indicator (1 = BCC, 2 = HCP)
    %
    % Outputs:
    %   schmidFactors    - Calculated raw Schmid factor matrix
    %   grainIDs1        - Flat array logging structural parent grain tracking IDs
    %   schmidFactors1   - Tracked mean Schmid factors mapped per rank entry
    %   slipSystems1     - Character string cell array displaying formatted Miller notation labels
    %   slipFamilies1    - Character string cell array mapping explicit HCP family descriptions
    %   traceAngles1     - Mapped raw angular offset measurements relative to horizon
    %   traceAngles1_fin - Mapped boundary-adjusted final trace projection angles
    %   q, qx, qy        - Graphic vector handles generated for plot management

    % Load crystal slip systems and variant counts for the active phase
    [ss, nList] = slip_systems(cs, m);

    % Establish an exclusion color palette avoiding green channels for colorblind clarity
    base_colors = [1 0 0; 0 0 1; 1 0 1; 1 1 0; 0 1 1; 1 0.5 0; 0.5 0 1];  
    colors = base_colors(mod(0:length(ss)-1, size(base_colors,1))+1, :);    

    % Initialize storage boundaries for nested loop performance optimization
    N = length(grains_phase);
    M = length(ss);
    row_idx = 1;
    grainIDs1      = NaN(N * M, 1);
    schmidFactors1 = NaN(N * M, 1);
    slipSystems1   = cell(N * M, 1);
    slipFamilies1  = cell(N * M, 1);
    traceAngles1   = NaN(N * M, 1);
    traceAngles1_fin   = NaN(N * M, 1);

    % Normalize macroscopic vector coordinate components to unified scale length
    stressDirection = normalize(loaddir);

    % -------------------------------------------------------------------------
    % GRAIN MATRIX LOOPING
    % -------------------------------------------------------------------------
    for g = 1:N
        % Isolate specific measurement pixels belonging to current index point
        grainEBSD = ebsd_phase(grains_phase(g));
        if isempty(grainEBSD)
            continue;
        end
        [~, grain_indices] = ismember(grainEBSD.id, ebsd_phase.id);
        orientations = grainEBSD.orientations;

        % Compute Schmid factors using local pixel orientation matrices
        schmidFactors = computeSchmidFactors(orientations, stressDirection, ss);
        [~, ss_active_phase] = max(schmidFactors, [], 2);

        % Set local operational variable mapping matching pixel elements
        ss_in_grain = ss_active_phase;

        % Calculate the arithmetic mean Schmid factor across all internal pixels
        meanSF_perSS = mean(schmidFactors, 1);  

        % Sort arrays in descending order to isolate premium mechanical load configurations
        [~, modes_per_rank] = sort(meanSF_perSS, 'descend');  

        % ---------------------------------------------------------------------
        % ITERATE SLIP SYSTEM VARIANTS WITHIN ACTIVE GRAIN
        % ---------------------------------------------------------------------
        for j = 1:M
            ss_idx = modes_per_rank(j);
            grainIDs1(row_idx)      = grains_phase(g).id;             
            schmidFactors1(row_idx) = mean(schmidFactors(:, ss_idx));             
            
            if m == 1 
                % -------------------------------------------------------------
                % PHASE 1: BETA (BCC) - 3-INDEX FORMATTING
                % -------------------------------------------------------------
                hkl = Miller(ss(ss_idx).n, cs{m}, 'plane');
                uvw = Miller(ss(ss_idx).b, cs{m}, 'direction');
                
                hkl_vals = round([hkl.h, hkl.k, hkl.l], 3);
                uvw_vals = round([uvw.h, uvw.k, uvw.l], 3); 
            
                ss_label = sprintf('(%d%d%d)[%d%d%d]', ...
                    hkl_vals(1), hkl_vals(2), hkl_vals(3), ...
                    uvw_vals(1), uvw_vals(2), uvw_vals(3));
                slipSystems1{row_idx} = ss_label;
                slipFamilies1{row_idx}  = '-';
                
                % Classify index group positioning to resolve color assignment
                if ss_idx <= nList(1)
                    color = colors(1, :);  
                elseif ss_idx <= nList(1)+nList(2)
                    color = colors(nList(1)+1, :); 
                else
                    color = colors(nList(1)+nList(2)+1, :); 
                end
            else 
                % -------------------------------------------------------------
                % PHASE 2: ALPHA (HCP) - 5-FAMILY ROBUST LOGIC
                % -------------------------------------------------------------
                limits = cumsum(nList);
                
                % Map specific system entry indices directly to structural labels
                if ss_idx <= limits(1)
                    ss_family = 'Prismatic';
                    cIdx = 1; 
                elseif ss_idx <= limits(2)
                    ss_family = 'Basal';
                    cIdx = limits(1) + 1;
                elseif ss_idx <= limits(3)
                    ss_family = 'Pyramidal <a>';
                    cIdx = limits(2) + 1;
                elseif ss_idx <= limits(4)
                    ss_family = 'Pyramidal I <c+a>';
                    cIdx = limits(3) + 1;
                else
                    ss_family = 'Pyramidal II <c+a>';
                    cIdx = limits(4) + 1; 
                end
                
                color = colors(cIdx, :);

                hkl = Miller(ss(ss_idx).n, cs{m}, 'plane');
                uvw = Miller(ss(ss_idx).b, cs{m}, 'direction');
                
                hkl_vals = round([hkl.h, hkl.k, hkl.i, hkl.l], 3);
                uvw_vals = round([uvw.U, uvw.V, uvw.T, uvw.W], 3); 
            
                ss_label = sprintf('(%d%d%d%d)[%d%d%d%d]', ...
                    hkl_vals(1), hkl_vals(2), hkl_vals(3), hkl_vals(4), ...
                    uvw_vals(1), uvw_vals(2), uvw_vals(3), uvw_vals(4));
                
                slipFamilies1{row_idx} = ss_family;
            end

            slipSystems1{row_idx} = ss_label;
            ori_mean = mean(orientations);
            slipPlane = ori_mean * ss(ss_idx).n;
        
            if isnan(slipPlane.x) || isnan(slipPlane.y) || isnan(slipPlane.z)
                continue;
            end
        
            % Compute spatial line intersection angles relative to horizontal frame
            [u, v, angle_deg, angle_deg_fin] = angle_with_horizontal(slipPlane);
            traceAngles1(row_idx) = angle_deg;
            traceAngles1_fin(row_idx) = angle_deg_fin;

            % Capture centroid point coordinates for graphic placement mapping
            c = grains_phase(g).centroid;

            % -----------------------------------------------------------------
            % GRAPHICAL OBJECT INJECTION: PLOT AND SLIP TRACE LINES
            % -----------------------------------------------------------------
            if j == 1
                % Top-ranked maximum system gets drawn as a bold solid vector
                q = quiver(c.x, c.y, u, v, 'Color', color, 'LineWidth', 1.5, 'MaxHeadSize', 1.5);
            else
                % Secondary candidate paths are drawn as dashed lines without headers
                q = quiver(c.x, c.y, u, v, 'Color', color, 'LineWidth', 1.5, 'LineStyle', '--' ,'MaxHeadSize', 0);
            end
            
            % Overlay base reference axes layers
            qx = quiver(c.x, c.y, 1, 0, 'Color', 'k', 'LineWidth', 0.5, 'MaxHeadSize', 0);
            qy = quiver(c.x, c.y, 0, 1, 'Color', 'k', 'LineWidth', 0.5, 'MaxHeadSize', 0);
            
            if j == 1
                % Append crisp, readable identification banners above top system origins
                text(c.x, c.y+0.5, num2str(grainIDs1(row_idx)), ...
                    'Color', color, 'FontSize', 10, 'FontWeight', 'bold', ...
                    'HorizontalAlignment', 'center','BackgroundColor', 'white', 'EdgeColor', 'black', 'Margin', 1);
            end
            row_idx = row_idx + 1;
        end
    end

    % -------------------------------------------------------------------------
    % EXPLICIT INTERFACE LEGEND STRUCTURING
    % -------------------------------------------------------------------------
    if m == 1
        legendLabels = ["{110}", "{112}", "{123}"];
        limits = cumsum(nList);
        legendIndices = [1, limits(1)+1, limits(2)+1];
    else 
        legendLabels = ["Prismatic", "Basal", "Pyramidal <a>", "Pyramidal I <c+a>", "Pyramidal II <c+a>"];
        limits = cumsum(nList);
        legendIndices = [1, limits(1)+1, limits(2)+1, limits(3)+1, limits(4)+1];
    end

    % Render zero-length vector hooks to map color legends without altering map data
    num_legs = length(legendLabels); 
    legendHandles = gobjects(1, num_legs);

    for h = 1:num_legs
        idx = legendIndices(h);
        if idx > size(colors, 1)
            idx = size(colors, 1); 
        end
        
        legendHandles(h) = quiver(0, 0, 0, 0, 'Color', colors(idx, :), 'LineWidth', 1.5, 'MaxHeadSize', 0);
    end

    % Finalize display box properties and layout constraints
    legend(legendHandles, legendLabels, ...
        'TextColor', 'black', ...
        'Location', 'northeastoutside', ...
        'FontSize', 10, ...
        'FontWeight', 'bold');
end