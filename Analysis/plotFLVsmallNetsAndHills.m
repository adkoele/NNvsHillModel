function plotFLVsmallNetsAndHills(A1, A2, x1, x2, mdlNames)

    datasets = {A1, A2};
    titles={'Force-Length Relationship', 'Force-Velocity Relationship'};
    xlabelS={'Length', 'Velocity'};
    xS={x1,x2};
    
    colors=lines(8);
    
    figure;
    tiledlayout(1,2);
    
    for gamma=1:2
        nexttile;
        hold on;
    
        for delta=1:length(mdlNames)
            plot(xS{gamma}, datasets{gamma}(:,delta), 'Color', colors(delta,:), 'LineWidth',1.5);
        end
    
        title(titles{gamma});
        xlabel(xlabelS{gamma});
        ylabel('Force');
        legend(mdlNames, 'Location', 'best');
    
    
    end
end