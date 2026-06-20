% =========================================================================
% MASTER SCRIPT: AUTOMATED MULTI-PHASE EBSD & SLIP TRACE ANALYSIS PIPELINE
% =========================================================================
%
% Purpose:
%   This top-level script automates the processing of electron backscatter
%   diffraction (.ctf) files for dual-phase Zirconium alloys (Zr-2.5% Nb).
%   It handles data calibration, phase separation (Alpha-HCP and Beta-BCC),
%   Schmid factor map evaluations, trace angle logging, and generates 
%   aggregated statistics on active crystallographic slip modes.

% Initialize the MTEX Toolbox environment
startup_mtex
clc; clear; close all;

% Define explicit system workspace directories for data sourcing and exports
inputDir = 'YOUR_INPUT_PATH_HERE';
outputBaseDir = 'YOUR_OUTPUT_PATH_HERE';

% Index all available .ctf datasets located within the input directory
ctfFiles = dir(fullfile(inputDir, '*.ctf'));

% Define fundamental crystal symmetries using precise lattice parameters
cS_hcp = crystalSymmetry('6/mmm', [3.26 3.26 5.13], 'mineral', 'Zr Nb Alpha HCP', 'X||a*', 'Y||b', 'Z||c');
cS_bcc = crystalSymmetry('m-3m', [3.555 3.555 3.555], 'mineral', 'Zr Nb Beta BCC', 'X||a*', 'Y||b', 'Z||c');
cs = {cS_bcc, cS_hcp}; % Combined crystal symmetry cell array

% Configure tracking metrics for the Alpha (HCP) phase slip families
labels_alpha = {'Prismatic', 'Basal', 'Pyramidal <a>', 'Pyramidal I <c+a>', 'Pyramidal II <c+a>'};
combinedSSCounts_alpha = zeros(1, numel(labels_alpha));
totalGrains_alpha = 0;

% Configure tracking metrics for the Beta (BCC) phase slip families
labels_beta = {'{110}<111>', '{112}<111>', '{123}<111>'};
combinedSSCounts_beta = zeros(1, numel(labels_beta));
totalGrains_beta = 0;

% =========================================================================
% LOOP THROUGH DETECTED DATASETS
% =========================================================================
for k = 1:length(ctfFiles)
    
    % Establish file paths and split names for organizing subfolders
    ctfPath = fullfile(ctfFiles(k).folder, ctfFiles(k).name);
    [~, sampleName, ~] = fileparts(ctfFiles(k).name);

    % Assign external stress/loading direction 
    loaddir = vector3d.X; % Modify as per the loading direction utilised
    z = vector3d.Z; % Out-of-plane reference vector

    % Generate and export the raw uncleaned crystal orientation map overlay image
    EBSD_crystal_orientation_image(ctfPath, cs, outputBaseDir, sampleName, z);
    
    % Execute full data cleaning, artifact filtering, and grain boundary smoothing
    [ebsd, grains] = EBSD_processing(ctfPath, cs);
    
    % Define 3D unit cell presentation dimensions and configurations for both phases
    for l = 2:3 
        cs0 = ebsd(ebsd.phaseId == l).CS;
        if l == 2
            cs_1 = crystalShape.cube(cs0);
            ipfKey_1 = ipfColorKey(cs0);
            ipfKey_1.inversePoleFigureDirection = z;
        else 
            cs_2 = crystalShape.hex(cs0);
            ipfKey_2 = ipfColorKey(cs0);
            ipfKey_2.inversePoleFigureDirection = z;
        end
    end
    
    % --- COMPOSITE MULTI-PHASE EBSD INTERACTION PLOT ---
    fig_ebsd = figure; 
    plot(ebsd(ebsd.phaseId == 2), ebsd(ebsd.phaseId == 2).orientations, ipfKey_1);
    hold on
    plot(ebsd(ebsd.phaseId == 3), ebsd(ebsd.phaseId == 3).orientations, ipfKey_2);
    hold on
    plot(grains.boundary, 'lineColor', [1 0 0], 'linewidth', 2); % Plot raw grain edges in red
    hold on
    plot(grains(grains.phaseId == 2), 0.4*cs_1, 'linewidth', 2, 'colored') % Superimpose cubic shapes
    hold on
    plot(grains(grains.phaseId == 3), 0.4*cs_2, 'linewidth', 2, 'colored') % Superimpose hexagonal shapes
    hold off
    
    set(gca, 'YDir', 'reverse'); % Match spatial image matrix coordinates
    legend off
    title('EBSD Map');

    % Create unique sample subdirectories and export the composite EBSD figure map
    outputFolder = fullfile(outputBaseDir, sampleName);
    if ~exist(outputFolder, 'dir'), mkdir(outputFolder); end
    set(fig_ebsd, 'Units', 'normalized', 'WindowState', 'maximized');
    drawnow;
    fileName = sprintf('%s_EBSD_Map.png', sampleName);
    exportgraphics(fig_ebsd, fullfile(outputFolder, fileName));
    close(fig_ebsd);

    % =========================================================================
    % PHASE SPECIFIC SUB-ANALYST LOOPS (m=1: Beta-BCC, m=2: Alpha-HCP)
    % =========================================================================
    for m = 1:2
        phaseLabel = ["Beta", "Alpha"]; 
        
        % Segment structural sub-properties targeted to current phase
        fig_ipdf = figure;
        ebsd_phase = ebsd(ebsd.phaseId == m+1);
        grains_phase = grains(grains.phaseId == m+1);
        
        % Plot Inverse Pole Figures relative to the surface plane normal vector
        plotIPDF(ebsd_phase.orientations, z);
        if m == 1
            title('Inverse Pole Figure - Beta Phase');
        else 
            title('Inverse Pole Figure - Alpha Phase');
        end

        % Maximize graphics canvas window and export Inverse Pole Figure
        set(fig_ipdf, 'Units', 'normalized', 'WindowState', 'maximized');
        drawnow;
        fileName = sprintf('%s_%s_IPDF.png', sampleName, phaseLabel(m));
        exportgraphics(fig_ipdf, fullfile(outputFolder, fileName));
        close(fig_ipdf);
        
        % --- GENERATE INVERSE POLE FIGURE COLOR SCHEME LEGEND KEYS ---
        fig_ipfkey = figure;
        if m == 1
            plot(ipfKey_1);
            title('IPF Key - Beta Phase (Z Direction)');
        else
            plot(ipfKey_2);
            title('IPF Key - Alpha Phase (Z Direction)');
        end

        set(fig_ipfkey, 'Units', 'normalized', 'WindowState', 'maximized');
        drawnow;
        fileName = sprintf('%s_%s_IPF_Key.png', sampleName, phaseLabel(m));
        exportgraphics(fig_ipfkey, fullfile(outputFolder, fileName));
        close(fig_ipfkey);
        
        % Extract map coloring data arrays corresponding to true spatial orientations
        if m == 1
            ipfColor = ipfKey_1.orientation2color(ebsd_phase.orientations);
        else
            ipfColor = ipfKey_2.orientation2color(ebsd_phase.orientations);
        end
        
        % --- CALCULATE AND OVERLAY THEORETICAL SURFACE SLIP TRACES ---
        fig_trace = figure; 
        plot(ebsd_phase, ipfColor);
        hold on;
        plot(grains.boundary, 'lineColor', 'k', 'lineWidth', 1.5); % Traced in black

        % Call out to internal function to project slip line options for every grain
        [schmidfactors, grainIDs, schmidFactors, slipSystems, slipFamilies, traceAngles, traceAngles_fin, q, qx, qy] = slip_trace(ebsd_phase, grains_phase, loaddir, cs, m);

        set(gca, 'YDir', 'reverse');
        title(sprintf('Slip Traces - %s Phase', phaseLabel(m))); 
        hold off;

        set(fig_trace, 'Units', 'normalized', 'WindowState', 'maximized');
        drawnow;
        fileName = sprintf('%s_%s_SlipTraces.png', sampleName, phaseLabel(m));
        exportgraphics(fig_trace, fullfile(outputFolder, fileName));
        close(fig_trace);

        % --- GENERATE EMBEDDED HIGH-VISIBILITY GRAIN ID VERIFICATION MAPS ---
        fig_grain_ids = figure;
        plot(ebsd_phase, ipfColor); 
        hold on;
        plot(grains_phase.boundary, 'lineColor', 'k', 'lineWidth', 1.5); 

        % Overlay text box identifiers precisely on grain geometric centers
        for i = 1:length(grains_phase)
            text(grains_phase(i).centroid.x, grains_phase(i).centroid.y, ...
                num2str(grains_phase(i).id), ...
                'FontSize', 8, 'FontWeight', 'bold', 'color', 'k', ...
                'HorizontalAlignment', 'center', 'BackgroundColor', 'w', 'EdgeColor', 'k', 'Margin', 1);
        end

        set(gca, 'YDir', 'reverse');
        title(sprintf('Grain ID Map - %s Phase', phaseLabel(m)));
        hold off;

        set(fig_grain_ids, 'Units', 'normalized', 'WindowState', 'maximized');
        drawnow;
        fileName = sprintf('%s_%s_GrainID_Map.png', sampleName, phaseLabel(m));
        exportgraphics(fig_grain_ids, fullfile(outputFolder, fileName));
        close(fig_grain_ids);

        % --- PACK AND SORT DATA MATRICES INTO EXCEL DATASHEETS ---
        valid = ~isnan(grainIDs); % Cull out unindexed data rows
        
        T = table(grainIDs(valid), schmidFactors(valid), slipSystems(valid), ...
                slipFamilies(valid), traceAngles(valid), traceAngles_fin(valid), ...
            'VariableNames', {'GrainID', 'SchmidFactor', 'SlipSystem', 'SlipFamily', 'TraceAngle_deg', 'TraceAngle_deg_fin'});

        outputFolder = fullfile(outputBaseDir, sampleName);
        if ~exist(outputFolder, 'dir'), mkdir(outputFolder); end
        
        tablePath = fullfile(outputFolder, sprintf('%s_%s_Grain_Slip_Data.csv', sampleName, phaseLabel(m)));
        writetable(T, tablePath);

        % Find and isolate unique grain labels to calculate max ranking states
        uniqueGrains = unique(T.GrainID);
        isMax = false(height(T), 1);

        % Flag the row containing the maximum Schmid factor for each grain
        for i = 1:length(uniqueGrains)
            thisGrain = uniqueGrains(i);
            rows = find(T.GrainID == thisGrain);

            if ~isempty(rows)
                [~, relIdx] = max(T.SchmidFactor(rows));
                isMax(rows(relIdx)) = true; % Toggle boolean flag high
            end
        end

        % Append logical column and write updated spreadsheet to disk
        T.IsMaxSchmid = isMax;
        writetable(T, tablePath);

        % --- HISTOGRAM PREPARATION AND RE-DISTRIBUTION ---
        [ss, ~] = slip_systems(cs, m);
        full_sf_map = computeSchmidFactors(ebsd_phase.orientations, loaddir, ss);
        ssCounts = slipsystemdist(full_sf_map, ebsd_phase, grains_phase, cs, m);

        % Store counts into globally aggregated variables for master graphing
        if m == 1
            combinedSSCounts_beta = combinedSSCounts_beta + ssCounts;
            totalGrains_beta = totalGrains_beta + length(grains_phase);
        else 
            combinedSSCounts_alpha = combinedSSCounts_alpha + ssCounts;
            totalGrains_alpha = totalGrains_alpha + length(grains_phase);
        end
    end
end
close all;

% =========================================================================
% FINAL DATA AGGREGATION & BAR GRAPH GENERATION
% =========================================================================

% Piece together individual arrays separated by an empty visual channel spacer
combinedCounts = [combinedSSCounts_alpha, 0, combinedSSCounts_beta]; 
xtickLabels = [labels_alpha, {''}, labels_beta]; 

figure;
bar(combinedCounts, 'FaceColor', [0.2 0.5 0.7]);
xticks(1:length(combinedCounts));
xticklabels(xtickLabels);
ylabel('Number of Grains with Active Slip System');
title('Aggregated Slip System Activation - \alpha & \beta Phases');
xtickangle(45);
grid on;

% Add floating analytical description boxes defining data boundaries
annotationText_alpha = sprintf('\\alpha - Total Grains: %d', totalGrains_alpha);
annotationText_beta = sprintf('\\beta - Total Grains: %d', totalGrains_beta);

annotation('textbox', [0.15 0.85 0.3 0.07], ...
    'String', annotationText_alpha, 'FontSize', 12, 'FontWeight', 'bold', ...
    'EdgeColor', 'black', 'BackgroundColor', 'white');

annotation('textbox', [0.55 0.85 0.3 0.07], ...
    'String', annotationText_beta, 'FontSize', 12, 'FontWeight', 'bold', ...
    'EdgeColor', 'black', 'BackgroundColor', 'white');

% Drop a vertical dashed separator line between the Alpha and Beta sections
xline(numel(labels_alpha) + 0.5, '--k', 'LineWidth', 1.2);

% Render graphics map canvas and export final summary document
set(gcf, 'Units', 'normalized', 'WindowState', 'maximized');
drawnow;
exportgraphics(gcf, fullfile(outputFolder, 'Combined_SlipSystem_vs_Grains_AlphaBeta.png'));

% --- FINAL CLEANUP AND RELEASE ---
fclose('all'); 
close all; 
disp('Processing Complete. All files released.');