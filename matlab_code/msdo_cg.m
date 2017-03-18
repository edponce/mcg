%% Multiple search directions with orthogonalization conjugate gradient
function [x,flag,relres,kiter,flop,arelres] = msdo_cg(A,b,x0,s,tol,kmax)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Authors:
%   Eduardo Ponce
%   Hartwig Anzt
%   Gregory D. Peterson
%
% The University of Tennessee, Knoxville
% October 2016
%
% Algorithm 1 MSDO-CG (page 14)
% Laura Grigori, Sophie Moufawad, Frederic Nataf. Enlarged Krylov Subspace Conjugate
% Gradient Methods for Reducing Communication. (Research Report) RR-8597,
% 2014. <hal-01065985>
%
% MSDO-CG Multiple Search Directions with Orthogonalization Conjugate Gradient method
%
% x = msdo_cg(A,b) return the approximate solution x for the system of linear equations A * x = b.
% The n-by-n coefficient matrix A must be symmetric positive definite matrix.
% The right-hand side column vector b must have length n.
%
% x = msdo_cg(A,b,x0) specifies the initial approximation for x using column vector
% x0 of size n. The default setting is a zero column vector.
%
% x = msdo_cg(A,b,x0,s) specifies the number of domain partitions. 
% The default setting is 1.
%
% x = msdo_cg(A,b,x0,s,tol) specifies the stopping tolerance of the method.
% The default setting is 1e-8.
%
% x = msdo_cg(A,b,x0,s,tol,kmax) specifies the maximum number of iterations
% of the method. The default setting is min(2 * n,1000).
%
% [x,flag] = msdo_cg(A,b) return information flag:
%   flag = 0: required tolerance satisfied
%   flag = 1: no convergence to the required tolerance within maximum
%             number of iterations
%   flag = 2: break down due to one of the iteration parameters becoming
%             zero/infinite
%   flag = 3: invalid domain partitioning
%
% [x,flag,relres] = msdo_cg(A,b) return the relative residual norm:
%   relres = ||b - A * x|| / ||b||.
%
% [x,flag,relres,kiter] = msdo_cg(A,b) return the number of iterations.
%
% [x,flag,relres,kiter,flops] = msdo_cg(A,b) return the number of
% floating-point operations required by the method.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Internal variables
n = size(A,1);  % dimension of A

% Output initializations
x = zeros(n,1);  % solution 
flag = 0;  % result flag
kiter = 0;  % iterations
flop = 0;  % FP cost

% NOTE: only for debugging
arelres = zeros(kmax,1);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Input A, nxn symmetric positive definite matrix

% Input b, nx1 right-hand side

% Input x0, nx1 initial guess or iterate 

% Input eps, stopping tolerance

% Input kmax, maximum number of iterations

% Output xk, approximate solution of A * x = b


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Cost = n - 1
% ||b||2
bnrm2 = norm(b);
if (bnrm2 == 0.0)
    bnrm2 = 1.0;
end

% 1. Cost = 2nnz + 2n - 1
% 1. r0 = b - A * x0
r = b - A * x0;

% 1. rho = (||r0||2)^2
relres = norm(r) / bnrm2;
arelres(1) = relres;
if (relres < tol)
    flag = 0;
    return;
end

% 1. k = 1
kiter = 1;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Domain partitions
if n < s
    flag = 3;
    fprintf('Error:\t%s\n', 'too many domains for matrix');
    return;
end
idxdoms = domain_partition(n,s);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Work variables
P = zeros(n,s);  % domain search directions
W = zeros(n,s);
T = zeros(n,s);  % temporary domain search directions
alpha = zeros(s,1);
beta = zeros(s,1);
Pis = zeros(n,s,20);
Wis = zeros(n,s,20);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 2. Cost = 2nnz - (t - 1)n
% 2. P1 = T(r0)
T = Tprojection(r, n, idxdoms, s);
P = T;

% 2. W1 = A * P1
W = A * P;
    

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 3. Cost = (6n - 1)(t - 1)(t/2) + (4n + 1)t
% 3. A-orthonormalize P1 against each others
[P,W] = Aorth_mgs_others(P, W, n, s);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 4. Cost = 2n
% 4. while (sqrt(rho) > eps * ||b||2 and k < kmax)
while (relres > tol && kiter < kmax)
    
    % Store Pi's and Wi's for A-orthonormalization
    Pis(1:n,1:s,kiter) = P;
    Wis(1:n,1:s,kiter) = W;

    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 5. Cost = (2n - 1)t
% 5. alpha = (Pk^t * Wk)^-1 * (Pk^t * r) = Pk^t * r
    alpha = P' * r;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 6. Cost = (2t - 1)n + n
% 6. x = x + Pk * alpha
    x = x + P * alpha;

    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 7. Cost = (2t - 1)n + n
% 7. r = r - Wk * alpha
    r = r - W * alpha;

    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 8. Cost = 2n - 1
% 8. rho = ||r||2^2 
    relres = norm(r) / bnrm2;
    arelres(kiter+1) = relres;

    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 9. Cost = (2n - 1)t
% 9. beta = -(Pk^t * Wk)^-1 * (Wk^t * r) = -Wk^t * r
    beta = -W' * r;
    

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 10. Cost = 2nt
% 10. Pk+1 = T(r) + Pk * diag(beta)
    T = Tprojection(r, n, idxdoms, s);
    P = T + P * diag(beta);

    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 11. Cost = 2nnz - (t - 1)n + 2nt
% 11. Wk+1 = A * T(rk) + Wk * diag(beta)
    W = A * T + W * diag(beta);
    

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 12. Cost = (6n - 1)kt^2 + (4n + 1)t
% 12. A-orthonormalize Pk+1 against all Pi's for i <= k
    %[P,W] = Aorth_mgs_previous(P, W, Pis, Wis, n, s, kiter);

    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 13. Cost = (6n - 1)(t - 1)(t/2) + (4n + 1)t
% 13. A-orthonormalize Pk+1 against each others
    [P,W] = Aorth_mgs_others(P, W, n, s);
    
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 14. Cost = 1
% 14. k = k+1
    kiter = kiter + 1;

% 15. end while
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Free storage
clear r

% Cost = 4nnz + 5n + k(11nt + 2nnz) = O(nnzk + ntk)
nnzA = nnz(A);
flop = 4 * nnzA + 5 * n + kiter * (11 * n * s + 2 * nnzA);
flop = flop + ((6 * n - 1) * (s - 1) * (s / 2) + (4 * n + 1) * s) * kiter;
flop = flop + ((6 * n - 1) * kiter * s^2 + (4 * n + 1) * s) * kiter;

% Check for breakdown
if isnan(relres)
    flag = 2;
end
end


% Compute indices for partitions of domains, search directions
function idxdoms = domain_partition(n,s)
    szdoms1= round(n / s);  % size of first/intermediate domains
    szdoms2 = n - ((s - 1) * szdoms1);  % size of last domain
    %szdoms2 = mod(n,szdoms1);  % size of last domain

    % Set domain indices
    idxdoms = zeros(s+1,1);  % array of domain indices
    idxdoms(1) = 1;  % first index
    for i = 2:s+1
        if i < s+1
            idxdoms(i) = idxdoms(i-1) + szdoms1;  % intermediate indices
        else
            idxdoms(i) = idxdoms(i-1) + szdoms2;  % last index
        end
    end
end


% Projection of vector on domain matrix (search directions)
function T = Tprojection(r, n, idxdoms, s)
    T = zeros(n,s);
    for i = 1:s
        ip1 = idxdoms(i);
        ip2 = idxdoms(i+1) - 1;
        T(ip1:ip2,i) = r(ip1:ip2); 
    end
end
