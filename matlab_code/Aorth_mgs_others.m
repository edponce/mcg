% Algorithm 13 (grigori p.38)
% A-orthonormalization of search directions Pk against each others
% Modified Gram-Schmidt
function [P, W] = Aorth_mgs_others(P, W, n, s)
    for i = 1:s
        for j = 1:i-1
            P(1:n,i) = P(1:n,i) - (W(1:n,j)' * P(1:n,i)) * P(1:n,j);
            W(1:n,i) = W(1:n,i) - (W(1:n,j)' * P(1:n,i)) * W(1:n,j);
        end
        pap = W(1:n,i)' * P(1:n,i);
%         if pap < 0
%             pap = -pap;
%         end
%         if pap < 1e-8
%             pap = 1.0;
%         end
        P(1:n,i) = P(1:n,i) / sqrt(pap);
        W(1:n,i) = W(1:n,i) / sqrt(pap);
    end
end
