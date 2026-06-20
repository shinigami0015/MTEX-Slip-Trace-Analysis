function cos_theta = cosTheta(U1,V1,W1,U2,V2,W2,a,c)
    % COSTHETA Calculates the cosine of the angle between two vectors in a Hexagonal system.
    %
    % Purpose:
    %   Unlike cubic systems where a simple vector dot product works, a hexagonal 
    %   (HCP) lattice is anisotropic (the 'a' and 'c' lattice parameters are unequal). 
    %   This function implements the generalized crystallographic metric tensor equation 
    %   to calculate the true geometric angle ($\theta$) between two directional vectors 
    %   using Miller-Bravais coordinate components [U, V, W] and the lattice constants.
    %
    % Inputs:
    %   U1, V1, W1 - Coordinates of the first vector in the hexagonal system
    %   U2, V2, W2 - Coordinates of the second vector in the hexagonal system
    %   a, c       - Lattice parameters (constants) of the hexagonal crystal structure
    %
    % Outputs:
    %   cos_theta  - The calculated cosine value ($\cos\theta$) of the angle between the two vectors

    % 1. Calculate the Numerator
    % This represents the modified dot product accounting for the 120-degree basal plane angle
    num = a^2 * ((U1*U2 + V1*V2) - 0.5*(U1*V2 + V1*U2)) ...
          + c^2 * W1*W2;

    % 2. Calculate the Denominator
    % This computes the magnitudes (lengths) of both vectors scaled by the lattice dimensions
    den = sqrt(a^2*(U1^2 - U1*V1 + V1^2) + c^2*W1^2) ...
          * sqrt(a^2*(U2^2 - U2*V2 + V2^2) + c^2*W2^2);

    % 3. Compute the Cosine of the Angle
    cos_theta = num / den;
end