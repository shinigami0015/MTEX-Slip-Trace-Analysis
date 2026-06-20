function EBSD_crystal_orientation_image (ctfPath, cs, outputBaseDir, sampleName, z)
    % EBSD_CRYSTAL_ORIENTATION_IMAGE Loads EBSD data, adjusts Euler angles, and saves maps.
    %
    % Purpose:
    %   This function imports raw electron backscatter diffraction (.ctf) data, 
    %   separates the material into its constituent phases (BCC and HCP), performs 
    %   Euler angle adjustments to align coordinates, generates spatial grain maps 
    %   overlaid with 3D crystal unit cell shapes, and exports high-resolution figures.
    %
    % Inputs:
    %   ctfPath       - Absolute system path to the input .ctf dataset
    %   cs            - Cell array containing crystal symmetry objects {BCC, HCP}
    %   outputBaseDir - Base folder path where processed folders will be saved
    %   sampleName    - Character string of the specific specimen being processed
    %   z             - Inverse Pole Figure (IPF) reference direction vector (typically vector3d.Z)

    % 1. Load EBSD Data 
    % Automatically converts Euler orientations to match spatial coordinates upon loading
    ebsd = EBSD.load(ctfPath, cs,'convertEuler2SpatialReferenceFrame');

    % Separate data by Phase ID (Phase 3 is Alpha-HCP, Phase 2 is Beta-BCC)
    ebsd_alpha = ebsd(ebsd.phaseId == 3);
    ebsd_beta = ebsd(ebsd.phaseId == 2);

    % 2. Extract Individual Euler Angles
    % Brakes orientations into Bunge-convention angles: phi1 (rotation 1), Phi (tilt), phi2 (rotation 2)
    [phi1_alpha, Phi_alpha, phi2_alpha] = Euler(ebsd_alpha.orientations);
    [phi1_beta, Phi_beta, phi2_beta] = Euler(ebsd_beta.orientations);

    % 3. Coordinate Alignment / Euler Modifications
    % Apply custom offset rotations to correct alignments between EBSD camera coordinates and physical samples
    phi1_alpha = mod(phi1_alpha + 90*degree, 360*degree);     % Add 90° offset to phi1 for Alpha
    phi2_alpha = mod(phi2_alpha - 30*degree, 360*degree);     % Subtract 30° offset from phi2 for Alpha

    phi1_beta = mod(phi1_beta + 90*degree, 360*degree);       % Add 90° offset to phi1 for Beta

    % Reconstruct modified orientations back into the EBSD master objects
    ebsd_alpha.orientations = orientation.byEuler(phi1_alpha, Phi_alpha, phi2_alpha, ebsd_alpha.CS);
    ebsd_beta.orientations = orientation.byEuler(phi1_beta, Phi_beta, phi2_beta, ebsd_beta.CS);

    % 4. Grain Reconstruction
    % Group neighboring pixels into grains based on a standard 10-degree misorientation threshold
    [grains, ebsd.grainId] = calcGrains(ebsd, 'angle', 10*degree);
    grains_alpha = grains(grains.phaseId == 3); 
    grains_beta = grains(grains.phaseId == 2); 

    grains_processed = grains;
    ebsd_processed = ebsd;

    % 5. Setup Inverse Pole Figure (IPF) Keys and 3D Crystal Overlays
    % Loop through Phase 2 (BCC) and Phase 3 (HCP) to prepare coloring keys and 3D shapes
    for l = 2:3 
        cs0 = ebsd(ebsd.phaseId == l).CS;
        
        if l == 2
            % Phase 2: Create a 3D Cube representation for Cubic symmetry
            cs_1 = crystalShape.cube(cs0);
            ipfKey_1 = ipfColorKey(cs0);
            ipfKey_1.inversePoleFigureDirection = z;
        else 
            % Phase 3: Create a 3D Hexagonal prism representation for Hexagonal symmetry
            cs_2 = crystalShape.hex(cs0);
            ipfKey_2 = ipfColorKey(cs0);
            ipfKey_2.inversePoleFigureDirection = z;
        end
    end

    % 6. Plotting the Microstructure Map
    figure; 
    
    % Plot Beta-BCC orientation map
    plot(ebsd(ebsd.phaseId == 2), ebsd(ebsd.phaseId == 2).orientations, ipfKey_1);
    hold on
    
    % Plot Alpha-HCP orientation map
    plot(ebsd(ebsd.phaseId == 3), ebsd(ebsd.phaseId == 3).orientations, ipfKey_2);
    hold on
    
    % Overlay 3D micro-cube indicators inside the BCC grain centers (scaled down to 40% size)
    plot(grains(grains.phaseId == 2), 0.4*cs_1, 'linewidth', 2, 'colored')
    hold on
    
    % Overlay 3D micro-hexagons inside the HCP grain centers (scaled down to 40% size)
    plot(grains(grains.phaseId == 3), 0.4*cs_2, 'linewidth', 2, 'colored')
    hold on
    
    % Trace solid black lines outlining all structural grain boundaries
    plot(grains.boundary, 'lineColor', 'k', 'lineWidth', 1.5);
    hold off
    
    % Final map cosmetic presentation configurations
    legend off
    title('EBSD Map');

    % 7. Auto-Save Figure to Disk
    % Locate the figure handle we just generated
    figHandles = findall(groot, 'Type', 'figure');

    % Check if destination directory exists, if not, construct it automatically
    outputFolder = fullfile(outputBaseDir, sampleName);
    if ~exist(outputFolder, 'dir')
        mkdir(outputFolder);
    end

    fig = figHandles;
    name = sprintf('%s_Figure_ebsd_orientations.png', sampleName);
    filePath = fullfile(outputFolder, name);

    try
        % Maximize the application window before exporting to guarantee crisp layout ratios
        figure(fig);  
        set(fig, 'Units', 'normalized', 'WindowState', 'maximized');
        drawnow;

        % Save figure to file path as a high-density image layout
        exportgraphics(fig, filePath);
    catch ME
        warning('Could not save figure %d: %s', i, ME.message);
    end

    % Clear active figure monitors to free up memory before moving to the next iteration
    close all;
end