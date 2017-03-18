function [P, J, Av, nz] = full2crs(A, baseIdx)
    
    m = size(A,1);          % Get matrix rows
    [nzi nzj nzv] = find(A);   % Get coordinates and values of non-zeros
    nz = length(nzv);       % Get number of non-zero elements
    P = zeros(m+1,1);       % Create row pointer array
    J = zeros(nz,1);        % Create column index array
    Av = zeros(nz,1);       % Create non-zero values array
 
    % Count elements on each row
    for i=1:nz
        P(nzi(i)+1) = P(nzi(i)+1)+1;
    end

    % Create row pointer by cumulative sum of elements
    P = cumsum(P);
    
    % Set column indices, values, and row pointers
    for i=1:nz
        Av(P(nzi(i))+1) = nzv(i);
        J(P(nzi(i))+1) = nzj(i);
        P(nzi(i)) = P(nzi(i))+1;
    end
    
    % Shift row pointer array
    for i=m:-1:1
        P(i+1) = P(i);
    end

    P(1) = 0;               % Set initial value of row pointer 
    if baseIdx > 0
        P = P + baseIdx;        % Add b-based indexing
    elseif baseIdx == 0
        J = J - 1;
    end
    
    clear m nzi nzj nzv;