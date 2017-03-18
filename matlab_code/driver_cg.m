%% Reset MATLAB environment
close all;
clear;
%clc;

%% Configurations

format long e;  % MATLAB numeric format

% User configurations
%infile = '../mm_matrix/Trefethen_2000.mtx';  % input matrix file
% infile = '../mm_matrix/crystm02.mtx';
%infile = '../mm_matrix/poisson2D/poisson2D.mtx';
%bfile = '../mm_matrix/poisson2D/poisson2D_b.mtx';  % input b vector

saveflag = 0;  % 0 = discard results, 1 = save results
methods = [1 3 4];  % 1 = cg, 2 = msdo_cg, 3 = lre_cg, 4 = msd_cg
plotflag = 1;  % 0 = no plots, 1 = result plots, 2 = result/intermediate plots

% Default settings
s = 1;  % number of search directions (domains)
tol = 1e-8;  % tolerance for convergence
kmax = 1000;  % maximum number of iterations


%% Initialization

% Set coefficient matrix, A
if exist('infile','var')
    [A, m, n, nz, repr, field, symm] = mmread(infile);
    fprintf('Input file:\t%s\n', infile);
    fprintf('Representation:\t%s\n', repr);
    fprintf('Field:\t%s\n', field);
    fprintf('Symmetry:\t%s\n', symm);
else
    A = gallery('poisson',100);
    [m, n] = size(A);
    nz = nnz(A);
end

% Validate matrix dimensions
if m ~= n
    fprintf('Error:\t%s\n', 'coefficient matrix is not symmetric');
end

% Validate matrix is not zero
if nz <= 0
    fprintf('Error:\t%s\n', 'coefficient matrix is zero');
end

% Create figure with nonzeros
if plotflag > 0
    spy(A);
    title('Nonzero diagram of A');
end

% Convert coefficient matrix to sparse format
% if ~strcmp(repr, 'array')
%     A = full(A);
% end

den = nz / (m * n);  % calculate matrix density

% Set initial LHS, x
%x0 = zeros(n, 1);
x0 = rand(n,1);

% Set initial RHS, b
if exist('bfile','var')
    [b, m, n, nz, repr, field, symm] = mmread(bfile);
else
    b = ones(m, 1);
end


% Testing framework variables
ntests = 4;
allrelres = zeros(kmax,ntests);  % store all residual vectors
allkiter = zeros(ntests,1);  % store all k iterations
testsflag = zeros(ntests,1);  % track tests run (for outputs)


%% Print outputs

% Print configurations and settings
fprintf('Tolerance:\t%g\n', tol);
fprintf('Max iterations:\t%d\n', kmax);
fprintf('A dims:\t%d x %d\n', size(A));
fprintf('b dims:\t%d x %d\n', size(b));
fprintf('x dims:\t%d x %d\n', size(x0));
fprintf('Matrix nonzeros:\t%d\n', nz);
fprintf('Matrix density:\t%g\n', den * 100);
fprintf('Matrix condition:\t%g\n', condest(A));
fprintf('\n');


%% CG

% Run method and measure runtime
tid = 1;
if any(methods == tid)
    testsflag(tid) = 1;
    fprintf('%s\n','CG Conjugate Gradient method');
    tic;
    [xcg, relres, kiter, flag, arelres] = cg(A, x0, b, kmax, tol);
    runtime = toc;
    
    % Calculate exact accuracy
    exactres = norm(b - A * xcg);

    % Print results
    fprintf('Flag:\t%d\n', flag);
    fprintf('Exact residual:\t%g\n', exactres);
    fprintf('Relative residual:\t%g\n', relres);
    fprintf('Iterations:\t%d\n', kiter);
    fprintf('Runtime:\t%f\n', runtime);
    fprintf('\n');
    allrelres(:,tid) = arelres;
    allkiter(tid) = kiter;
end


%% MSDO-CG

% Run method and measure runtime
tid = 2;
if any(methods == tid)
    testsflag(tid) = 1;
    fprintf('%s\n','MSDO-CG Multiple Search Directions with Orthogonalization Conjugate Gradient method');
    tic;
    [xmsdo, flag, relres, kiter, flop, arelres] = msdo_cg(A, b, x0, s, tol, kmax);
    runtime = toc;

    % Calculate exact accuracy
    exactres = norm(b - A * xmsdo);

    % Calculate flops
    flops = flop / runtime;

    % Print results
    fprintf('Flag:\t%d\n', flag);
    fprintf('Domains:\t%d\n', s);
    fprintf('Exact residual:\t%g\n', exactres);
    fprintf('Relative residual:\t%g\n', relres);
    fprintf('Iterations:\t%d\n', kiter);
    fprintf('Runtime:\t%f\n', runtime);
    fprintf('Flops:\t%g\n', flops);
    fprintf('\n');
    allrelres(:,tid) = arelres;
    allkiter(tid) = kiter;
end


%% LRE-CG

% Run method and measure runtime
tid = 3;
if any(methods == tid)
    testsflag(tid) = 1;
    fprintf('%s\n','LRE-CG Long Recurrence Enlarged Conjugate Gradient method');
    tic;
    [xlre, flag, relres, kiter, flop, arelres] = lre_cg(A, b, x0, s, tol, kmax);
    runtime = toc;

    % Calculate exact accuracy
    exactres = norm(b - A * xlre);

    % Calculate flops
    flops = flop / runtime;

    % Print results
    fprintf('Flag:\t%d\n', flag);
    fprintf('Domains:\t%d\n', s);
    fprintf('Exact residual:\t%g\n', exactres);
    fprintf('Relative residual:\t%g\n', relres);
    fprintf('Iterations:\t%d\n', kiter);
    fprintf('Runtime:\t%f\n', runtime);
    fprintf('Flops:\t%g\n', flops);
    fprintf('\n');
    allrelres(:,tid) = arelres;
    allkiter(tid) = kiter;
end


%% MSD-CG

% Run method and measure runtime
tid = 4;
if any(methods == tid)
    testsflag(tid) = 1;
    fprintf('%s\n','MSD-CG Multiple Search Directions Conjugate Gradient method');
    tic;
    [xmsd, flag, relres, kiter, flop, arelres] = msd_cg(A, b, x0, s, tol, kmax);
    runtime = toc;

    % Calculate exact accuracy
    exactres = norm(b - A * xmsd);

    % Calculate flops
    flops = flop / runtime;

    % Print results
    fprintf('Flag:\t%d\n', flag);
    fprintf('Domains:\t%d\n', s);
    fprintf('Exact residual:\t%g\n', exactres);
    fprintf('Relative residual:\t%g\n', relres);
    fprintf('Iterations:\t%d\n', kiter);
    fprintf('Runtime:\t%f\n', runtime);
    fprintf('Flops:\t%g\n', flops);
    fprintf('\n');
    allrelres(:,tid) = arelres;
    allkiter(tid) = kiter;
end


% NOTE: plot for debugging
if plotflag > 0
    plot_residuals(allkiter, allrelres);
end
