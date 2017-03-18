%% Long recurrence enlarged conjugate gradient
function [x,flag,relres,kiter,flop,arelres] = lre_cg(A,b,x0,s,tol,kmax)
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
% Algorithm 2 LRE-CG (page 16)
% Laura Grigori, Sophie Moufawad, Frederic Nataf. Enlarged Krylov Subspace Conjugate
% Gradient Methods for Reducing Communication. (Research Report) RR-8597,
% 2014. <hal-01065985>
%
% LRE-CG Long Recurrence Enlarged Conjugate Gradient method
%
% x = lre_cg(A,b) return the approximate solution x for the system of linear equations A * x = b.
% The nxn coefficient matrix A must be symmetric positive definite matrix.
% The right-hand side column vector b must have length n.
%
% x = lre_cg(A,b,x0) specifies the initial approximation for x using column vector
% x0 of size n. The default setting is to use a zero column vector.
%
% x = lre_cg(A,b,x0,s) specifies the number of domain partitions. 
% The default setting is 1.
%
% x = lre_cg(A,b,x0,s,tol) specifies the stopping tolerance of the method.
% The default setting is to use 1e-8.
%
% x = lre_cg(A,b,x0,s,tol,kmax) specifies the maximum number of iterations
% of the method. The default setting is min(2xn,1000).
%
% [x,flag] = lre_cg(A,b) return information flag:
%   flag = 0: required tolerance satisfied
%   flag = 1: no convergence to the required tolerance within maximum
%             number of iterations
%   flag = 2: break down due to one of the iteration parameters becoming
%             zero/infinite
%
% [x,flag,relres] = lre_cg(A,b) return the relative residual norm:
%   relres = ||b - A * x|| / ||b||.
%
% [x,flag,relres,kiter] = lre_cg(A,b) return the number of iterations.
%
% [x,flag,relres,kiter,flops] = lre_cg(A,b) return the number of
% floating-point operations required by the method.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Internal variables
n = size(A, 1);
        
% Output initializations
x = zeros(n, 1);  % solution 
flag = 0;  % flag result
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
Q = zeros(n,s);  % changes in size
W = zeros(n,s);
T = zeros(n,s);  % temporary domain search directions


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 2. Cost = 3n
% 2. W = T(r0)
T = Tprojection(r, n, idxdoms, s);
W = T;

% 2. Q = normalize(W)
W = normcols(W, n, s);  % normalize columns
Q = W;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 3. Cost = 2n
% 3. while (sqrt(rho) > eps * ||b||2 and k < kmax)
while (relres > tol && kiter < kmax)
    

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 4. Cost = (2nnz - n) * t * k + (2n - 1)t^2 * k^2
% 4. G = Q' * A * Q
    ks = kiter * s;
    G = zeros(ks,ks);  % changes in size
    G = Q' * A * Q;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 5. Cost = (2n - 1)t * k + Solve_alpha(t * k)
% 5. alpha = G^-1 * (Q' * r)
    alpha = inv(G) * (Q' * r);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 6. Cost = 2tkn
% 6. x = x + Q * alpha
    x = x + Q * alpha;
    

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 7. Cost = 2tkn
% 7. r = r - A * Q * alpha
    r = r - A * Q * alpha;
    

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 8. Cost = 2n - 1
% 8. rho = ||r||2^2 
    relres = norm(r) / bnrm2;
    arelres(kiter+1) = relres;

    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 9. Cost = (2nnz - n)t
% 9. W = A * W
    W = A * W;
    

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 12. Cost = ?
% 12. Orthonormalize W against Q
    W = orth_mgs_others(W, Q, n, s, kiter);
    

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 13. Cost = ?
% 13. Orthonormalize W
% 13. Q = [Q W]
    W = normcols(W, n, s);
    
    Qk = Q;
    Q = zeros(n,(kiter+1)*s);
    Q = [Qk W];
    
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 14. Cost = 1
% 14. k = k+1
    kiter = kiter + 1;

% 15. end while
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


% Normalize columns of matrix
function W = normcols(W, n, s)
    for i = 1:s
        fac = norm(W(1:n,i));
        W(1:n,i) = W(1:n,i) / fac;
    end
end


% Orthonormalization of vectors W against the vectors of Q
% Modified Gram-Schmidt
function W = orth_mgs_others(W, Q, n, s, k)
    ks = k * s;
    for i = 1:s
        for j = 1:ks
            W(1:n,i) = W(1:n,i) - (Q(1:n,j)' * W(1:n,i)) * Q(1:n,j);
        end
        pap = W(1:n,i)' * W(1:n,i);
%         if pap < 0
%             pap = -pap;
%         elseif pap == 0
%             pap = 1.0;
%         end
        W(1:n,i) = W(1:n,i) / sqrt(pap);
    end
end
