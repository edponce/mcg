% NOTE: debugging, plot residuals
function plot_residuals(allkiter, allrelres)
    figure;
    hold on;
    [n, nt] = size(allrelres);
    mymarkers = {'k+-','r*','bo','m.-'};
    mytests = {'cg','msdo','lre','msd'};
    lgnd = [];
    for i = 1:nt
        k = allkiter(i);
        if (k > 0)
            x = 1:k;
            y = allrelres(1:k,i);
            semilogy(x, y, mymarkers{i});
            %plot(x, y, mymarkers(i));
            lgnd = [lgnd {mytests{i}}];
        end
    end
    title('CG relative residuals');
    legend(lgnd);
    hold off;
end