function [ss, nList] = slip_systems(cs, m)
    % SLIP_SYSTEMS Generates symmetrical slip systems and counts for BCC and HCP phases.
    %
    % Purpose:
    %   This function serves as the central crystallographic database for the 
    %   deformation pipeline. Depending on the targeted phase, it defines the 
    %   primary slip planes (n) and slip directions (b), applies antipodal 
    %   symmetrization to account for forward/backward slip equivalence, 
    %   and records the total number of variants available within each slip family.
    %
    % Inputs:
    %   cs    - Cell array containing crystal symmetry objects {BCC, HCP}
    %   m     - Current active phase index indicator (1 = Beta-BCC, 2 = Alpha-HCP)
    %
    % Outputs:
    %   ss    - Vertical column vector containing all symmetric slipSystem objects
    %   nList - Row vector storing the count of individual variants per slip family

    if m == 1
        % -----------------------------------------------------------------
        % PHASE 1: BETA (BCC) - CUBIC SLIP FAMILIES
        % -----------------------------------------------------------------
        
        % Define {110}<111> slip system family with antipodal symmetry
        ss110 = slipSystem(Miller({1, 1, 0}, cs{m}), Miller({1, -1, 1}, cs{m}), cs{m});
        ss110 = ss110.symmetrise('antipodal');

        % Define {112}<111> slip system family with antipodal symmetry
        ss112 = slipSystem(Miller({1, 1, 2}, cs{m}), Miller({1, 1, -1}, cs{m}), cs{m});
        ss112 = ss112.symmetrise('antipodal');

        % Define {123}<111> slip system family with antipodal symmetry
        ss123 = slipSystem(Miller({1, 2, 3}, cs{m}), Miller({1, 1, -1}, cs{m}), cs{m});
        ss123 = ss123.symmetrise('antipodal');

        % Log the discrete variant counts for each cubic slip family
        nList = [length(ss110), length(ss112), length(ss123)];

        % Concatenate all active systems into a single vertical vector block
        ss = [ss110; ss112; ss123];
        
    else
        % -----------------------------------------------------------------
        % PHASE 2: ALPHA (HCP) - HEXAGONAL SLIP FAMILIES
        % -----------------------------------------------------------------
        
        % Extract native MTEX HCP slip configurations using antipodal symmetry
        ssPrism  = slipSystem.prismaticA(cs{m}).symmetrise('antipodal');
        ssBasal  = slipSystem.basal(cs{m}).symmetrise('antipodal');
        ssPyr    = slipSystem.pyramidalA(cs{m}).symmetrise('antipodal');
        ssPyrCA1 = slipSystem.pyramidalCA(cs{m}).symmetrise('antipodal');  
        ssPyr2   = slipSystem.pyramidal2CA(cs{m}).symmetrise('antipodal');

        % Force column vectors (semicolon stacking format)
        ss = [ssPrism; ssBasal; ssPyr; ssPyrCA1; ssPyr2];

        % Ensure nList variant counts match the vertical array assembly exactly
        nList = [length(ssPrism), length(ssBasal), length(ssPyr), length(ssPyrCA1), length(ssPyr2)];
    end
end