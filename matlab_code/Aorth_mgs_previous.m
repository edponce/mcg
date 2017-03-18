% Algorithm 12 (grigori p.38)
% A-orthonormalization of search directions Pk against previous Pi's
% Modified Gram-Schmidt
function [P, W] = Aorth_mgs_previous(P, W, Pis, Wis, n, s, k)
    for o = 1:s
        for i = 1:k
            Pi = Pis(1:n,1:s,i);
            Wi = Wis(1:n,1:s,i);
            for j = 1:s
                P(1:n,o) = P(1:n,o) - (Wi(1:n,j)' * P(1:n,o)) * Pi(1:n,j);
                W(1:n,o) = W(1:n,o) - (Wi(1:n,j)' * P(1:n,o)) * Wi(1:n,j);
            end
        end
        pap = W(1:n,o)' * P(1:n,o);
%         if pap < 0
%             pap = -pap;
%         end
%         if pap < 1e-8
%             pap = 1.0;
%         end
        P(1:n,o) = P(1:n,o) / sqrt(pap);
        W(1:n,o) = W(1:n,o) / sqrt(pap);
    end
end
