function plotSmallNNLosses()
    %% Check path
    filenamePath=mfilename('fullpath');
    filePath =[fileparts(filenamePath) filesep '..' filesep];
    addpath(genpath(filePath))

    model_params = {'network_r01_15_05_2025', ...
    'network_r12_15_05_2025'};


    for alpha=1:length(model_params)
        load(model_params{alpha});

        for beta=1:length(res)
            network(alpha).iterationNum{beta}=res(beta).Mdl.TrainingHistory.Iteration;
            network(alpha).trainingLoss{beta}=res(beta).Mdl.TrainingHistory.TrainingLoss;
            network(alpha).validationLoss{beta}=res(beta).Mdl.TrainingHistory.ValidationLoss;
        end
    end

    modelNames={'NN-r01', 'NN-r12'};
    colors={[0.6 0.8 1],[1 0.6 0.6]};
    figure;
    tiledlayout(1, length(model_params));

    for alpha=1:length(model_params)
        maxLenTrain = max(cellfun(@length, network(alpha).trainingLoss));
        maxLenValid = max(cellfun(@length, network(alpha).validationLoss));

        trainMat = zeros(length(model_params), maxLenTrain);
        valMat = zeros(length(model_params), maxLenTrain);

        for beta=1:length(res)
            t=network(alpha).trainingLoss{beta};
            v=network(alpha).validationLoss{beta};
            t=t(:)';
            v=v(:)';
            trainMat(beta,:)=[t, repmat(t(end),1,maxLenTrain-length(t))];
            valMat(beta,:)=[v, repmat(v(end),1,maxLenTrain-length(v))];
        end

        trainMean = mean(trainMat,1);
        trainStd=std(trainMat,0,1);
        valMean=mean(valMat,1);
        valStd=std(valMat,0,1);

        x=1:maxLenTrain;

        axmain = nexttile;
        hold on;

        for beta=1:length(res)
            lgdTr(beta)=plot(x, trainMat(beta,:), 'Color', [colors{1} 0.8], 'LineWidth',1);
            lgdVal(beta)=plot(x, valMat(beta,:), 'Color', [colors{2} 0.8], 'LineWidth',1);
        end

        lgdMeanTr(alpha)=plot(x, trainMean, 'Color', colors{1}*0.6, 'LineWidth',3);        
        lgdMeanVal(alpha)=plot(x, valMean, 'Color', colors{2}*0.6, 'LineWidth',3);
        title([modelNames{alpha} '- Loss Curves']);
        xlabel('Iteration');
        ylabel('Loss');
        legend([lgdTr(10), lgdVal(10), lgdMeanTr(alpha), lgdMeanVal(alpha)], ...
            {'Individual network - Training','Individual network - Validation', ...
            'Mean - Training','Mean - Validation'},'Location','northeast', 'Box','off');


        if alpha==1
            xlimZoom=[30 40];
        else
            xlimZoom=[55 70];
        end

        pos = axmain.Position;
        insetWidth=0.3*pos(3);
        insetHeight=0.3*pos(4);
        insetLeft=pos(1)+0.6*pos(3);
        insetBottom=pos(2)+0.55*pos(4);

        axInset=axes('Position', [insetLeft insetBottom insetWidth insetHeight]);
        hold on;

        idx = x>= xlimZoom(1) & x<= xlimZoom(2);

        plot(x(idx), trainMean(idx), 'Color', colors{1}*0.6, 'LineWidth',3);
        plot(x(idx), valMean(idx), 'Color', colors{2}*0.6, 'LineWidth',3);
        zoomLength = xlimZoom(2)-xlimZoom(1);
        title(['Zoomed view of last ' num2str(zoomLength) ' iterations']);

        axes(axmain);

        
    end

end