%% Store matrix to file
function storeCRS(fn, A, P, J, m, n, nz)

  fd = fopen(fn, 'wt');           % Create a text file to save matrix data.
  fprintf(fd, '%d\t%d\t%d\n', m, n, nz);
  for i = 1:max(m+1,nz)
    if i <= m+1 && i <= nz
        fprintf(fd, '%d\t%d\t%f\n', P(i), J(i), A(i));   % Write matrix data
    elseif i <= m+1
        fprintf(fd, '%d\n', P(i));            % Write matrix data
    else
        fprintf(fd, '\t%d\t%f\n', J(i), A(i));   % Write matrix data
    end
  end
  fclose(fd);

  clear fd