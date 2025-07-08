function plotBigNNFLVs()

    %% Check if force-length and force-velocity relationships are followed for best network
    nnFil='Pu1_all_all_5network';
    model_names = {'NN', 'Hill'};
    model_params={nnFil,'normal'};
    muscle_names = {'LG', 'DF'};
    folder_networksparams = '/path/to/code/'; %Change this to the folder where the networks and model parameter sets are stored

    
    act = 0.2:0.2:1;
    
    colors = {'#f3e79b', '#fac484', '#eb7f86', '#ce6693', '#5c53a5'};
    styles = {'-', '--'};
    width = 2;
    
    figure

    for i = 1:length(act)
        l_ce1 = 0.6:0.01:1.4;
        l_ce1 = l_ce1';
        for j = 1:length(model_names)
                flce(:,j,i) = getForce(folder_networksparams,l_ce1, zeros(size(l_ce1)), act(i)+zeros(size(l_ce1)), model_names{j}, model_params{j}, muscle_names{1});
        end
        
        v_ce1 = -10:0.1:10;
        v_ce1 = v_ce1';
        for j = 1:length(model_names)
                gvce(:,j,i) = getForce(folder_networksparams,ones(size(v_ce1)), v_ce1, act(i)+zeros(size(v_ce1)), model_names{j}, model_params{j}, muscle_names{1});
        end
    end


    maxFlce=zeros(1,length(act));
    lceAtMax=cell(1,length(act));

    for i=1:length(act)
        forceVec=flce(:,1,i);
        maxFlce(i)=max(forceVec);
        idx=find(forceVec==maxFlce(i));
        lceAtMax{i}=l_ce1(idx);
    end

    maxGlce=zeros(1,length(act));
    vceAtMax=cell(1,length(act));

    for i=1:length(act)
        forceVec2=gvce(:,1,i);
        maxGlce(i)=max(forceVec2);
        idx2=find(forceVec2==maxGlce(i));
        vceAtMax{i}=v_ce1(idx2);
    end




    
    subplot(1,2,1)
    hold on
    for i = 1:length(act)
        for j = 1:length(model_names)
            plot(l_ce1, flce(:,j,i), "Color", colors{i}, "LineStyle", styles{j}, "LineWidth", width)

            if strcmp(model_names{j}, 'NN')
                [~,idxMax]=max(flce(:,j,i));
                xMax=l_ce1(idxMax);
                yMax=flce(idxMax,j,i);
                plot(xMax, yMax, 'kd', 'MarkerFaceColor','k', 'MarkerSize',10, 'HandleVisibility','off');
            end
        end
    end
    xlabel('Normalized fibre length')
    ylabel('Normalized muscle force')
    title('Force-Length Relationship')

    
    mdls = {'NN', 'Hill'};

    subplot(1,2,2)
    hold on
    plotCounter = 1;
    legendEntries=cell(1,length(act)*length(model_names));
    for i = 1:length(act)
        for j = 1:length(model_names)
            plot(v_ce1, gvce(:,j,i), "Color", colors{i}, "LineStyle", styles{j}, "LineWidth", width)
            
            legendEntries{plotCounter} = [mdls{j}, ', a=' num2str(act(i))];
            plotCounter = plotCounter + 1;
            if strcmp(model_names{j}, 'NN')
                [~,idxMax]=max(gvce(:,j,i));
                xMax=v_ce1(idxMax);
                yMax=gvce(idxMax,j,i);
                plot(xMax, yMax, 'kd', 'MarkerFaceColor','k', 'MarkerSize',10, 'HandleVisibility','off');
            end
        end
    end

    hDummy=plot(NaN, NaN, 'kd', 'MarkerFaceColor','k', 'MarkerSize',10);
    legendEntries{end+1}='Max Force';
    legend([legendEntries, {'Max Force'}]);
    xlabel('Normalized fibre velocity')
    ylabel('Normalized muscle force')
    title('Force-Velocity Relationship')





    %% Plot FLVs for remaining networks

    nnFil = {'Bl3_all_all_4network', 'BL4_all_all_4network', 'Or3_all_all_2network'};
    model_names = {'NN', 'Hill'};
    
    act = 0.2:0.2:1;
    
    colors = {'#f3e79b', '#fac484', '#eb7f86', '#ce6693', '#5c53a5'};
    styles = {'-', '--'};
    width = 2;
    
    figure
    
    for f=1:3
        model_params={nnFil{f},'normal'};
        for i = 1:length(act)
            l_ce1 = 0.6:0.01:1.4;
            l_ce1 = l_ce1';
            for j = 1:length(model_names)
                    flce(:,j,i) = getForce(folder_networksparams,l_ce1, zeros(size(l_ce1)), act(i)+zeros(size(l_ce1)), model_names{j}, model_params{j}, muscle_names{1});
            end
            
            v_ce1 = -10:0.1:10;
            v_ce1 = v_ce1';
            for j = 1:length(model_names)
                    gvce(:,j,i) = getForce(folder_networksparams,ones(size(v_ce1)), v_ce1, act(i)+zeros(size(v_ce1)), model_names{j}, model_params{j}, muscle_names{1});
            end
        end
    
    

        maxFlce=zeros(1,length(act));
        lceAtMax=cell(1,length(act));
    
        for i=1:length(act)
            forceVec=flce(:,1,i);
            maxFlce(i)=max(forceVec);
            idx=find(forceVec==maxFlce(i));
            lceAtMax{i}=l_ce1(idx);
        end
    
        maxGlce=zeros(1,length(act));
        vceAtMax=cell(1,length(act));
    
        for i=1:length(act)
            forceVec2=gvce(:,1,i);
            maxGlce(i)=max(forceVec2);
            idx2=find(forceVec2==maxGlce(i));
            vceAtMax{i}=v_ce1(idx2);
        end


    
        
        subplot(2,3,f)
        hold on
        for i = 1:length(act)
            for j = 1:length(model_names)
                plot(l_ce1, flce(:,j,i), "Color", colors{i}, "LineStyle", styles{j}, "LineWidth", width)
                if strcmp(model_names{j}, 'NN')
                    [~,idxMax]=max(flce(:,j,i));
                    xMax=l_ce1(idxMax);
                    yMax=flce(idxMax,j,i);
                    plot(xMax, yMax, 'kd', 'MarkerFaceColor','k', 'MarkerSize',10, 'HandleVisibility','off');
                end
            end
        end
        xlabel('Normalized fibre length')
        ylabel('Normalized muscle force')
        title('Force-Length Relationship')
    
        
        mdls = {'NN', 'Hill'};
    
        subplot(2,3,f+3)
        hold on
        plotCounter = 1;
        legendEntries=cell(1,length(act)*length(model_names));
        for i = 1:length(act)
            for j = 1:length(model_names)
                plot(v_ce1, gvce(:,j,i), "Color", colors{i}, "LineStyle", styles{j}, "LineWidth", width)
                
                legendEntries{plotCounter} = [mdls{j}, ', a=' num2str(act(i))];
                plotCounter = plotCounter + 1;
                if strcmp(model_names{j}, 'NN')
                    [~,idxMax]=max(gvce(:,j,i));
                    xMax=v_ce1(idxMax);
                    yMax=gvce(idxMax,j,i);
                    plot(xMax, yMax, 'kd', 'MarkerFaceColor','k', 'MarkerSize',10, 'HandleVisibility','off');
                end
            end
        end
        hDummy=plot(NaN, NaN, 'kd', 'MarkerFaceColor','k', 'MarkerSize',10);
        legendEntries{end+1}='Max Force';
        legend([legendEntries, {'Max Force'}]);
        xlabel('Normalized fibre velocity')
        ylabel('Normalized muscle force')
        title('Force-Velocity Relationship')
    end

    text(0.02, 0.75, 'Force-Length', 'Units', 'normalized', ...
        'Rotation', 90, 'FontWeight', 'bold', 'FontSize', 18, ...
        'HorizontalAlignment','center');
     text(0.02, 0.25, 'Force-Velocity', 'Units', 'normalized', ...
        'Rotation', 90, 'FontWeight', 'bold', 'FontSize', 18, ...
        'HorizontalAlignment','center');

end