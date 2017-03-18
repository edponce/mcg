%% Store matrix to file
function storeCCS(fn, A, I, P, m, n, nz)

  fd = fopen(fn, 'wt');           % Create a text file to save matrix data.
  fprintf(fd, '%d\t%d\t%d\n', m, n, nz);
  for i = 1:max(m+1,nz)
    if i <= m+1 && i <= nz
        fprintf(fd, '%d\t%d\t%f\n', I(i), P(i), A(i));   % Write matrix data
    elseif i <= m+1
        fprintf(fd, '\t%d\n', P(i));            % Write matrix data
    else
        fprintf(fd, '%d\t\t%f\n', I(i), A(i));   % Write matrix data
    end
  end
  fclose(fd);

  clear fd