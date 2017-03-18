function [I, J, Av, nz] = full2aij(A, baseIdx)
   
    [I, J, Av] = find(A);       % Get coordinates and values of non-zeros
    nz = length(Av);            % Get number of non-zero elements
    
    % Convert to 0-based indexing if necessary
    if baseIdx == 0
        I = I - 1;
        J = J - 1;
    end