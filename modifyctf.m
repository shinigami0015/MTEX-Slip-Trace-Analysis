% =========================================================================
% SCRIPT: AUTOMATED CTF FILE PHASE RE-INDEXING & RE-ALIGNMENT TOOL
% =========================================================================
%
% Purpose:
%   This utility fixes an issue where the EBSD software incorrectly maps phase IDs 
%   (e.g., swapping Phase 1 and Phase 2 indices). It reads a raw text-based .ctf file,
%   finds the structural header where phase definitions are established, swaps the 
%   definition order, and then re-indexes the Phase ID column (Column 1) for every
%   spatial measurement coordinate before exporting a corrected file.

% Target absolute file path selection
filename = "D:\Abhinav Chandraker (Pls do not delete)\Zr alloy\Zr slip trace\All CTF\Longitudinal\160425 10% strain A2 Samp4 BARC S5.ctf";

% Read the full text payload into a vertical string matrix array
lines = readlines(filename);

% Locate the specific text line sequence marking the phase property setup boundaries
phaseLineIdx = find(contains(lines, 'Phases'));

% Extract the total integer count of structural mineral phases declared in the file
nPhases = sscanf(lines(phaseLineIdx), 'Phases %d');

% Locate the physical row locations of Phase 1 and Phase 2 definitions
phase1Line = phaseLineIdx + 1;
phase2Line = phaseLineIdx + 2;

% Swap the entire string definitions between the lines to correct structural indexing
temp = lines(phase1Line);
lines(phase1Line) = lines(phase2Line);
lines(phase2Line) = temp;

% Calculate the row index where individual pixel coordinate blocks start
dataStart = phaseLineIdx + nPhases + 1;

% -------------------------------------------------------------------------
% ITERATIVE PIXEL DATA PASS: PHASE VALUE RE-MAPPING
% -------------------------------------------------------------------------
for i = dataStart+1:length(lines)
    if strlength(lines(i)) < 5  % Skip empty space gaps and end-of-file artifacts
        continue;
    end
    
    % Split the current line string based on blank white-space delimiters
    tokens = split(lines(i));