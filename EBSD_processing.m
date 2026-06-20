function [ebsd_processed, grains_processed] = EBSD_processing(ctfPath, cs)
    % EBSD_PROCESSING Imports, calibrates, filters, and reconstructs EBSD data.
    %
    % Purpose:
    %   This function handles the core data-cleaning pipeline for raw data. It 
    %   loads the data, applies angular correction offsets, filters out non-physical 
    %   outlier grains (noise/artifacts) based on specific phase size thresholds, 
    %   smoothes pixel fluctuations using a spline filter, and computes final polished 
    %   grain boundary maps.
    %
    % Inputs:
    %   ctfPath          - Absolute system path to the input .ctf dataset
    %   cs               - Cell array containing crystal symmetry objects {BCC, HCP}
    %
    % Outputs:
    %   ebsd_processed   - Cleaned, smoothed MTEX EBSD dataset
    %   grains_processed - Polished MTEX grain object ready for downstream slip trace analysis

    % 1. Load EBSD Data 
    % Forces data alignment between Euler orientations and physical spatial coordinates upon loading
    ebsd = EBSD.load(ctfPath, cs, 'convertEuler2SpatialReferenceFrame');

    % Separate data data points by Phase ID (Phase 3 is Alpha-HCP, Phase 2 is Beta-BCC)
    ebsd_alpha = ebsd(ebsd.phaseId == 3);
    ebsd_beta = ebsd(ebsd.phaseId == 2);

    % 2. Extract Individual Euler Angles
    % Breaks orientations into standard Bunge-convention angles: phi1, Phi, phi2
    [phi1_alpha, Phi_alpha, phi2_alpha] = Euler(ebsd_alpha.orientations);
    [phi1_beta, Phi_beta, phi2_beta] = Euler(ebsd_beta.orientations);

    % 3. Coordinate Realignment / Calibration
    % Adjust orientations to correct angular offsets between scanning reference axes and sample geometry
    phi1_alpha = mod(phi1_alpha + 90*degree, 360*degree);     % Add 90° to phi1 for Alpha
    phi2_alpha = mod(phi2_alpha - 30*degree, 360*degree);     % Subtract 30° from phi2 for Alpha

    phi1_beta = mod(phi1_beta + 90*degree, 360*degree);       % Add 90° to phi1 for Beta

    % Reassign corrected orientations back to their respective phase structures
    ebsd_alpha.orientations = orientation.byEuler(phi1_alpha, Phi_alpha, phi2_alpha, ebsd_alpha.CS);
    ebsd_beta.orientations = orientation.byEuler(phi1_beta, Phi_beta, phi2_beta, ebsd_beta.CS);

    % 4. Initial Grain Finding for Data Purging
    % Run a fast grain reconstruction check at a 10-degree threshold to locate size-dependent noise
    [grains, ebsd.grainId] = calcGrains(ebsd, 'angle', 10*degree);
    grains_alpha = grains(grains.phaseId == 3); 
    grains_beta = grains(grains.phaseId == 2); 

    % 5. Noise and Artifact Removal
    % Delete unphysical, extremely small Alpha grains (less than 20 pixels) often caused by indexing noise
    ebsd(grains_alpha(grains_alpha.grainSize < 20)) = [];
    % Delete excessively large Beta grains (greater than 100 pixels) representing indexing errors or artifacts
    ebsd(grains_beta(grains_beta.grainSize > 100)) = [];

    % 6. Fine Grain Boundary Calculation
    % Re-calculate grains from the purged indexed dataset using a refined 2-degree threshold 
    % and a strict minimum requirement of 5 pixels per grain.
    [grains, ebsd('indexed').grainId] = calcGrains(ebsd('indexed'), 'alpha', 2, 'angle', 2*degree, 'minPixel', 5);
    
    % 7. Map Smoothing and Void Filling
    % Apply an MTEX spline filter to smooth indexed orientation data map-wide and fill in unindexed gaps
    ebsd = smooth(ebsd('indexed'), splineFilter, 'fill', grains);
    
    % Re-verify final grain shapes after spatial orientation smoothing has occurred
    [grains, ebsd('indexed').grainId] = calcGrains(ebsd('indexed'), 'alpha', 2, 'angle', 2*degree, 'minPixel', 5);

    % 8. Boundary Spline Curvature Polishing
    % Smooth out jagged pixel steps along the calculated grain boundaries for a natural morphology
    grains = smooth(grains, 4);

    % 9. Return Polished Outputs
    grains_processed = grains;
    ebsd_processed = ebsd;
end