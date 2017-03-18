%% Store matrix to file
function storeAIJ(fn, A, I, J, m, n, nz)

  fd = fopen(fn, 'wt');           % Create a text file to save matrix data.
  fprintf(fd, '%d\t%d\t%d\n', m, n, nz);
  for i = 1:nz
        fprintf(fd, '%d\t%d\t%f\n', I(i), J(i), A(i));   % Write matrix data
  end
  fclose(fd);

  clear fd