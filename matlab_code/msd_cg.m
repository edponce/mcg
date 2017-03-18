%% Multiple search directions conjugate gradient
function [x,flag,relres,kiter,flop,arelres] = msd_cg(A,b,x0,s,tol,kmax)
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
% Algorithm 2.2 MSD-CG (page 1137)
% Tongxiang Gu, Xingping Liu, Zeyao Mo, Xuebin Chi. Multiple search
% direction conjugate gradient method I: methods and their propositions.
% International Journal of Computer Mathematics. 2004.
%
% Algorithm 9 MSD-CG (page 36)
% Laura Grigori, Sophie Moufawad, Frederic Nataf. Enlarged Krylov Subspace Conjugate
% Gradient Methods for Reducing Communication. (Research Report) RR-8597,
% 2014. <hal-01065985>
%
% MSD-CG Multiple Search Directions Conjugate Gradient method
%
% x = msd_cg(A,b) return the approximate solution x for the system of linear equations A * x = b.
% The n-by-n coefficient matrix A must be symmetric positive definite matrix.
% The right-hand side column vector b must have length n.
%
% x = msd_cg(A,b,x0) specifies the initial approximation for x using column vector
% x0 of size n. The default setting is a zero column vector.
%
% x = msd_cg(A,b,x0,s) specifies the number of domain partitions. 
% The default setting is 1.
%
% x = msd_cg(A,b,x0,s,tol) specifies the stopping tolerance of the method.
% The default setting is 1e-8.
%
% x = msd_cg(A,b,x0,s,tol,kmax) specifies the maximum number of iterations
% of the method. The default setting is min(2 * n,1000).
%
% [x,flag] = msd_cg(A,b) return information flag:
%   flag = 0: required tolerance satisfied
%   flag = 1: no convergence to the required tolerance within maximum
%             number of iterations
%   flag = 2: break down due to one of the iteration parameters becoming
%             zero/infinite
%
% [x,flag,relres] = msd_cg(A,b) return the relative residual norm:
%   relres = ||b - A * x|| / ||b||.
%
% [x,flag,relres,kiter] = msd_cg(A,b) return the number of iterations.
%
% [x,flag,relres,kiter,flops] = msd_cg(A,b) return the number of
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
% ||b||2
bnrm2 = norm(b);
if (bnrm2 == 0.0)
    bnrm2 = 1.0;
end

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
T = zeros(n,s);  % temp domain search directions
alpha = zeros(s,1);
beta = zeros(s,1);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 2-3. P = T(r0)
T = Tprojection(r, n, idxdoms, s);
P = T;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 4. while (sqrt(rho) > eps * ||b||2 and k < kmax)
while (relres > tol && kiter < kmax)
    
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 5. alpha = (P^t * A * P)^-1 * (P^t * r)
    alpha = inv(P' * A * P) * (P' * r);
    

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 6. x = x + P * alpha
    x = x + P * alpha;

    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 7. r = r - A * P * alpha
    r = r - A * P * alpha;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 8. beta = -(P^t * A * P)^-1 * (P^t * A * r)
    beta = -inv(P' * A * P) * (P' * A * r);
    

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 9-10. P(:,i) = T(r) + beta(i) * P(:,i)
    T = Tprojection(r, n, idxdoms, s);
    P = T + P * diag(beta);

    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 11. rho = ||r||2^2 
    relres = norm(r) / bnrm2;
    arelres(kiter+1) = relres;
    
      
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 12. k = k+1
    kiter = kiter + 1;

% 13. end while
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Free storage
clear r

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
