function ssCounts = slipsystemdist(schmidfactors, ebsd_phase, grain_phase, cs, m)
    % SLIPSYSTEMDIST Categorizes and plots active slip family statistics across all grains.
    %
    % Purpose:
    %   This statistical function evaluates the dominant deformation mode for an entire
    %   phase map. It identifies the maximum active slip system index per EBSD pixel, 
    %   determines the statistical mode for each grain boundary group, classifies these 
    %   modes into their respective crystallographic slip families (3 for BCC, 5 for HCP), 
    %   and renders a polished histogram bar plot summarizing grain-by-grain activation.
    %
    % Inputs:
    %   schmidfactors - Precomputed global matrix tracking Schmid factors map-wide
    %   ebsd_phase    - Isolated phase subset of the processed MTEX EBSD dataset
    %   grain_phase   - Extracted structural grain handles for the target phase
    %   cs            - Cell array containing crystal symmetry objects {BCC, HCP}
    %   m             - Current active phase index indicator (1 = BCC, 2 = HCP)
    %
    % Outputs:
    %   ssCounts      - Flat row matrix capturing the final sum tally per slip family group

    % -------------------------------------------------------------------------
    % STEP 1: COMPUTE MAXIMUM SLIP FACTORS FOR ALL SCANNED ORIENTATIONS
    % -------------------------------------------------------------------------
    [~, nList] = slip_systems(cs, m);
    [~, ss_active_all] = max(schmidfactors, [], 2);  % Isolate peak index position for each EBSD point

    % -------------------------------------------------------------------------
    % STEP 2: ARBITRATE DOMINANT SLIP SYSTEM PER INDIVIDUAL GRAIN
    % -------------------------------------------------------------------------
    ss_modes = zeros(length(grain_phase),1);
    for g = 1:length(grain_phase)
        grain_ebsd = ebsd_phase(grain_phase(g));
        [~, grain_indices] = ismember(grain_ebsd.id, ebsd_phase.id);
        ss_in_grain = ss_active_all(grain_indices);
        ss_modes(g) = mode(ss_in_grain); % Calculate dominant variant trend via statistical mode
    end

    % -------------------------------------------------------------------------
    % STEP 3: CATEGORIZE DISCRETE SLIP SYSTEMS INTO CRYSTALLOGRAPHIC FAMILIES
    % -------------------------------------------------------------------------
    if m == 1
        % --- PHASE 1: BETA (BCC) CUBIC DEFORMATION FAMILIES ---
        labels = {'110', '112', '123'};
        ssCounts = zeros(1, 3);
        range110 = 1:nList(1);  % ss110
        range112 = (range110(end)+1):(range110(end)+nList(2));  % ss112
        range123 = (range112(end)+1):(range112(end)+nList(3));  % ss123

        % Map matched systems to cubic group arrays
        ssCounts(1) = sum(ismember(ss_modes, range110));
        ssCounts(2) = sum(ismember(ss_modes, range112));
        ssCounts(3) = sum(ismember(ss_modes, range123));
    else 
        % --- PHASE 2: ALPHA (HCP) HEXAGONAL DEFORMATION FAMILIES ---
        labels = {'Prismatic', 'Basal', 'Pyramidal <a>', 'Pyr I <c+a>', 'Pyr II <c+a>'};
        ssCounts = zeros(1, 5);

        % Establish consecutive indexing bounds using cumulative limits
        rangePrism = 1:nList(1);
        rangeBasal = (rangePrism(end)+1):(rangePrism(end)+nList(2));
        rangePyr = (rangeBasal(end)+1):(rangeBasal(end)+nList(3));
        rangePyrCA1 = (rangePyr(end)+1):(rangePyr(end)+nList(4));      
        rangePyr2   = (rangePyrCA1(end)+1):(rangePyrCA1(end)+nList(5)); 

        % Accumulate final frequency numbers for each hexagonal family
        ssCounts(1) = sum(ismember(ss_modes, rangePrism));
        ssCounts(2) = sum(ismember(ss_modes, rangeBasal));
        ssCounts(3) = sum(ismember(ss_modes, rangePyr));
        ssCounts(4) = sum(ismember(ss_modes, rangePyrCA1));
        ssCounts(5) = sum(ismember(ss_modes, rangePyr2));
    end

    % -------------------------------------------------------------------------
    % STEP 4: FIGURES AND FREQUENCY PLOT RENDERING
    % -------------------------------------------------------------------------
    figure;
    bar(ssCounts, 'FaceColor', [0.2 0.4 0.6]);
    xticks(1:length(labels));
    xticklabels(labels);
    ylabel('Number of Grains with Active Slip System');
    title('Slip Traces vs Slip System');
    grid on;

    % -------------------------------------------------------------------------
    % STEP 5: OVERLAY METADATA STATISTICS ANNOTATION BANNER
    % -------------------------------------------------------------------------
    totalGrains = length(grain_phase);
    annotationText = sprintf('Total Grains: %d', totalGrains);
    annotation('textbox', [0.15 0.75 0.2 0.1], ...
        'String', annotationText, ...
        'FontSize', 15, ...
        'FontWeight', 'bold', ...
        'EdgeColor', 'black', ...
        'HorizontalAlignment', 'center', ...
        'VerticalAlignment', 'middle', ...
        'BackgroundColor', 'white', ...
        'LineWidth', 1.2);

    hold off
end