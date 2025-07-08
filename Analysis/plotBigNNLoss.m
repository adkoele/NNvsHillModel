function plotBigNNLoss()
    %% Check path
    filenamePath=mfilename('fullpath');
    filePath =[fileparts(filenamePath) filesep '..' filesep];
    addpath(genpath(filePath))
    folder_networksparams = '/path/to/code/'; 


    %this filename pattern should match the filename of the best NN (lowest
    %error value). In this case it was the neural network type which
    %excluded bird 5 (Pu1) from training
    filenamePatterns={'Pu1*network.mat'};
    mdlFileArray={};
    
    for i=1:length(filenamePatterns)
        files=dir(fullfile(folder_networksparams,filenamePatterns{i}));
        fullNames=fullfile({files.folder}, {files.name});
        mdlFileArray = [mdlFileArray, fullNames];
    end
    
    mdlFileArray=mdlFileArray(:);


    for alpha=1:length(mdlFileArray)
        load(mdlFileArray{alpha});
        network(alpha).iterationNum=NN.Mdl.TrainingHistory.Iteration;
        network(alpha).trainingLoss=NN.Mdl.TrainingHistory.TrainingLoss;
        network(alpha).validationLoss=NN.Mdl.TrainingHistory.ValidationLoss;       
    end

    modelNames={'NN-b5'};
    colors={[0.6 0.8 1],[1 0.6 0.6]};
    figure;
    
    lossLengths=arrayfun(@(x) length(x.trainingLoss), network);
    maxLenTrain=max(lossLengths);


        trainMat = zeros(length(mdlFileArray), maxLenTrain);
        valMat = zeros(length(mdlFileArray), maxLenTrain);

        for beta=1:length(mdlFileArray)
            t=network(beta).trainingLoss;
            v=network(beta).validationLoss;
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

        
        hold on;

        for beta=1:length(mdlFileArray)
            lgdTr(beta)=plot(x, trainMat(beta,:), 'Color', [colors{1} 0.8], 'LineWidth',1);
            lgdVal(beta)=plot(x, valMat(beta,:), 'Color', [colors{2} 0.8], 'LineWidth',1);
        end

        lgdMeanTr=plot(x, trainMean, 'Color', colors{1}*0.6, 'LineWidth',3);        
        lgdMeanVal=plot(x, valMean, 'Color', colors{2}*0.6, 'LineWidth',3);
        title([modelNames '- Loss Curves']);
        xlabel('Iteration');
        ylabel('Mean Squared Error');
        legend([lgdTr(5), lgdVal(5), lgdMeanTr, lgdMeanVal], ...
            {'Individual network - Training','Individual network - Validation', ...
            'Mean - Training','Mean - Validation'},'Location','northeast', 'Box','off');

        %% Plot supporting figure of other 3 big network types
        filenamePatterns={'Bl3*network.mat', 'BL4*network.mat','Or3*network.mat'};
        modelLabels={'NN-b1', 'NN-b2', 'NN-b3'};

        figure;
        tiledlayout(1, length(filenamePatterns),'TileSpacing','compact','Padding','compact');

        for alpha=1:length(filenamePatterns)
            files=dir(fullfile(folder_networksparams, filenamePatterns{alpha}));
            fullNames=fullfile({files.folder},{files.name});
            mdlFileArray=fullNames(:);

            for beta=1:length(mdlFileArray)
                load(mdlFileArray{beta});
                network(beta).iterationNum=NN.Mdl.TrainingHistory.Iteration;
                network(beta).trainingLoss=NN.Mdl.TrainingHistory.TrainingLoss;
                network(beta).validationLoss=NN.Mdl.TrainingHistory.ValidationLoss;
            end

            lossLengths=arrayfun(@(x) length(x.trainingLoss), network);
            maxLenTrain=max(lossLengths);

            trainMat = zeros(length(mdlFileArray), maxLenTrain);
            valMat = zeros(length(mdlFileArray), maxLenTrain);
    
            for beta=1:length(mdlFileArray)
                t=network(beta).trainingLoss(:)';
                v=network(beta).validationLoss(:)';
                
                trainMat(beta,:)=[t, repmat(t(end),1,maxLenTrain-length(t))];
                valMat(beta,:)=[v, repmat(v(end),1,maxLenTrain-length(v))];
            end
    
            trainMean = mean(trainMat,1);
            valMean=mean(valMat,1);
            x=1:maxLenTrain;

            nexttile;
            hold on;

            hIndTrain = [];
            hIndVal = [];
            hMeanTrain = [];
            hMeanVal = [];
            for beta=1:length(mdlFileArray)
                h1=plot(x,trainMat(beta,:),'Color',[colors{1} 0.8],'LineWidth',1);
                h2=plot(x,valMat(beta,:),'Color',[colors{2} 0.8],'LineWidth',1);
                if beta==2
                    hIndTrain=h1;
                    hIndVal=h2;
                end

            end

            hMeanTrain=plot(x, trainMean,'Color',colors{1}*0.6,'LineWidth',2.5);
            hMeanVal=plot(x, valMean,'Color',colors{2}*0.6,'LineWidth',2.5);

            title([modelLabels{alpha} ' - Loss Curves']);

            
            ylabel('Mean Squared Error');
            box off;

            if alpha==2
                xlabel('Iteration');
                legend([hIndTrain, hIndVal, hMeanTrain, hMeanVal], {'Individual network - Training','Individual network - Validation', ...
            'Mean - Training','Mean - Validation'}, 'Location','best','Box','off');
            end
        end
  
end