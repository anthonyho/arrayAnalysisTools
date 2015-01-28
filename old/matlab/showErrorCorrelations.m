function showErrorCorrelations(data,errors)
% SHOWERRORCORRELATIONS Shows the correlation plot of a variables and its 
% error
%
% Anthony Ho, 10/28/2014


    %% Plotting correlation plot
    
    diagLine = linspace(min(data),max(data),100000);
    
    figure;
    hold on;
    
    semilogy(data,errors,'ob');
    makepretty(1,7);
    semilogy(diagLine,diagLine,'-r','LineWidth',1.5);
    
    hold off;
    
    set(gca,'yscale','log')

    
end