clear all
close all
clc

%% Add to path
%Please make sure that you use the correct \ or / for your operating system
addpath(genpath('/path/to/code/')) %Change this to the folder where you downloaded this code

%% Settings

bird_names = {'Bl3', 'BL4'};
muscle_names = {'LG', 'DF'};

%Write the location of the data and networks/model parameters here. Please
%make sure that you use the correct \ or / for your operating system
folder_data = 'path/to/data/'; %Change this to the folder where the guinnea fowl data is available
folder_networsparams = 'path/to/networks/and/params/files/'; %Change this to the folder where the networks and model parameter sets are stored

bird_data = readmatrix([folder_data 'MuscleMorphologyData']);
warning('The columns of the xlsx file are hard coded, please make sure that they match your version of MuscleMorphologyData.xlsx in getMuscleParameters.m')



model_params = {'optresult_05062025_Bl3d2_r01_baselineHillModel.mat', ...
    'optresult_05062025_Bl3d2_r12_baselineHillModel.mat', ...
    'network_r01_15_05_2025', ...
    'network_r12_15_05_2025'};


model_names = {'Hill', 'Hill', 'NN', 'NN' };
modelNamesWTrial={'Hill-r01', 'Hill-r12', 'NN-r01', 'NN-r12'};



%% Plot small NN losses
plotSmallNNLosses();


%% Check FL an FV relationships
%===========================================
%% Check if force-length and force-velocity relationships are followed
model_paramsFLV = {'optresult_05062025_Bl3d2_r01_baselineHillModel.mat', ...
    'optresult_05062025_Bl3d2_r12_baselineHillModel.mat', ...
    'network_r01_15_05_2025', ...
    'network_r12_15_05_2025'};
model_namesFLV = {'Hill', 'Hill', 'NN', 'NN'};
actFLV = 1;
for i = 1:length(actFLV)
    l_ce1 = 0.6:0.001:1.4;
    l_ce1 = l_ce1';
    for j = 1:length(model_namesFLV)
            flce(:,j,i) = ...
            getForce(folder_networsparams, ...
            l_ce1, zeros(size(l_ce1)), ...
            actFLV(i)+zeros(size(l_ce1)), ...
            model_namesFLV{j}, model_paramsFLV{j}, ...
            muscle_names{1});
    end
    
    v_ce1 = -10:0.01:10;
    v_ce1 = v_ce1';
    for j = 1:length(model_namesFLV)
            gvce(:,j,i) = ...
            getForce(folder_networsparams, ...
            ones(size(v_ce1)), v_ce1, ...
            actFLV(i)+zeros(size(v_ce1)), ...
            model_namesFLV{j}, model_paramsFLV{j}, ...
            muscle_names{1});
    end
end

modelNamesFLVplots={'Hill-r01', 'Hill-r12', 'NN-r01', 'NN-r12'};
plotFLVsmallNetsAndHills(flce, gvce, l_ce1, v_ce1, modelNamesFLVplots);

colors = {'#f3e79b', '#fac484', '#eb7f86', '#ce6693', '#5c53a5'};
styles = {'-', '--'};
width = 2;

figure
subplot(1,2,1)
hold on
for i = 1:length(actFLV)
    for j = 1:length(model_names)
        plot(l_ce1, flce(:,j,i), "Color", colors{i}, "LineStyle", styles{j}, "LineWidth", width)
    end
end
xlabel('Normalized fibre length')
ylabel('Normalized muscle force')
title('Force-Length Relationship')

plotCounter = 1;
mdls = {'NN', 'Hill'};

subplot(1,2,2)
hold on
for i = 1:length(actFLV)
    for j = 1:length(model_names)
        plot(v_ce1, gvce(:,j,i), "Color", colors{i}, "LineStyle", styles{j}, "LineWidth", width)
        txt=[mdls{j}, ', a=' num2str(actFLV(i))];
        legendEntries{plotCounter} = txt;
        plotCounter = plotCounter + 1;
    end
end
legend(legendEntries)
xlabel('Normalized fibre velocity')
title('Force-Velocity Relationship')

%===========================================




ModelNameCol = [];
BirdCol = [];
MuscleCol = [];

%Metric Columns
avgRMSE = [];
stdRMSE = [];
avgCorr = [];
medianR2 = [];
minR2 = [];
maxR2 = [];
meanMissedPcnt = [];
meanExtraPcnt = [];
minMissedPcnt = [];
maxMissedPcnt = [];
minExtraPcnt = [];
maxExtraPcnt = [];
rmseRise = [];
rmseFall = [];

R2counterToPlot=0;

for k = 1:length(bird_names)
    trial_names = dir([folder_data filesep bird_names{k}]);
    % Find some indices for average calculations and plotting
    if strcmpi(bird_names{k}, 'bl3')
        clear trial_name
        trial_names = dir([folder_data filesep bird_names{k}]); 
        for iT = 1:length(trial_names)
            trial_name{iT} = trial_names(iT).name;
        end
        r01_ind = find(~contains(trial_name, 'r01') == 0);
        r10_ind = find(~contains(trial_name, 'r10') == 0);
        r12_ind = find(~contains(trial_name, 'r12') == 0);
    elseif strcmpi(bird_names{k}, 'bl4')
        clear trial_name
        trial_names = dir([folder_data filesep bird_names{k}]); 
        for iT = 1:length(trial_names)
            trial_name{iT} = trial_names(iT).name;
        end
        r09_ind = find(~contains(trial_name, 'r09') == 0);
    end

    results(k).birdID=bird_names{k};

    timeAtMetricGenStart = datetime('now');
    disp(['Start of metric generation: ', char(timeAtMetricGenStart)]);

    for l = 1:length(muscle_names) 
        results(k).muscles(l).muscleName=muscle_names{l};
        for j = 1:length(model_names)
            results(k).muscles(l).models(j).modelName=modelNamesWTrial{j};

            allRMSE = [];
            allR2 = [];
            allCorr = [];
            allMissed = [];
            allExtra = [];
            allMeasured = [];
            allPredicted = [];
            allRiseRMSE = [];
            allFallRMSE = [];

            R2counterToPlot=R2counterToPlot+1;
            
            for i = 3:length(trial_names)
                if ~strcmp(trial_names(i).name(end-2:end), 'mat')
                    continue
                end 

                if strcmp(bird_names{k}, 'Bl3') && strcmp(muscle_names{l}, 'LG') && contains(model_params{j}, 'r01') && i==r01_ind
                    disp(['Skipping predictions for trial 01 for bird ', num2str(k), ' muscle ', num2str(muscle_names{l}), ' and network ', num2str(model_params{j})]);
                    continue
                end

                if strcmp(bird_names{k}, 'Bl3') && strcmp(muscle_names{l}, 'LG') && contains(model_params{j}, 'r12') && i==r12_ind
                    disp(['Skipping predictions for trial 12 for bird ', num2str(k), ' muscle ', num2str(muscle_names{l}), ' and network ', num2str(model_params{j})]);                    
                    continue
                end

                if strcmp(bird_names{k}, 'Bl3') && strcmp(muscle_names{l}, 'DF') && contains(model_params{j}, 'r01') && i==r01_ind && contains(model_params{j}, 'widthOptOnDF') 
                    disp(['Skipping predictions for trial 01 for bird ', num2str(k), ' muscle ', num2str(muscle_names{l}), ' and network ', num2str(model_params{j})]);                    
                    continue
                end
                if strcmp(bird_names{k}, 'Bl3') && strcmp(muscle_names{l}, 'DF') && contains(model_params{j}, 'r12') && i==r12_ind && contains(model_params{j}, 'widthOptOnDF')
                    disp(['Skipping predictions for trial 12 for bird ', num2str(k), ' muscle ', num2str(muscle_names{l}), ' and network ', num2str(model_params{j})]);                    
                    continue
                end

            
                %% Get forces, calculate RMSEs and correlations
                [time(i,k,l).t, data(i,k,l).l_ce, data(i,k,l).v_ce, data(i,k,l).EMG, force(i,1,k,l).meas] = loadDataFile(bird_names{k}, muscle_names{l}, trial_names(i).name,folder_data); 
                                          
                results(k).muscles(l).models(j).trials(i).measuredForce=force(i,1,k,l).meas;
                results(k).muscles(l).models(j).trials(i).time=time(i,k,l).t;

                if startsWith(model_params{j}, 'network')
                    

                    nnFile=load([folder_networsparams model_params{j}]);
                    
                    
                    N=length(nnFile.res);
                    
                    for iInd=1:length(nnFile.res)

                        [force(i,j,k,l).est,f_max(i,j,k,l)] = getForce(folder_networsparams, data(i,k,l).l_ce, data(i,k,l).v_ce, data(i,k,l).EMG, model_names{j}, model_params{j}, muscle_names{l}, bird_names{k}, bird_data, iInd);

                        results(k).muscles(l).models(j).trials(i).nnVariant(iInd).nnVariantID=iInd;
                        results(k).muscles(l).models(j).trials(i).nnVariant(iInd).estimatedForce=force(i,j,k,l).est;
                        results(k).muscles(l).models(j).trials(i).nnVariant(iInd).fMax=f_max(i,j,k,l);

                        [allRMSE, allCorr, allR2, allMeasured, allPredicted, allMissed, allExtra, allRiseRMSE, allFallRMSE] = getMetrics(i, j, k, l, force, f_max, iInd, muscle_names, time, allRMSE, allCorr, allR2, allMeasured, allPredicted, allMissed, allExtra, allRiseRMSE, allFallRMSE);

                    end
                    
                else
                    
                    [force(i,j,k,l).est,f_max(i,j,k,l)] = getForce(folder_networsparams, data(i,k,l).l_ce, data(i,k,l).v_ce, data(i,k,l).EMG, model_names{j}, model_params{j}, muscle_names{l}, bird_names{k}, bird_data);                    

                    results(k).muscles(l).models(j).trials(i).modelVariant(1).modelVariantID=1;
                    results(k).muscles(l).models(j).trials(i).modelVariant(1).estimatedForce=force(i,j,k,l).est;
                    results(k).muscles(l).models(j).trials(i).modelVariant(1).fMax=f_max(i,j,k,l);

                    [allRMSE, allCorr, allR2, allMeasured, allPredicted, allMissed, allExtra, allRiseRMSE, allFallRMSE] = getMetrics(i, j, k, l, force, f_max, i, muscle_names, time, allRMSE, allCorr, allR2, allMeasured, allPredicted, allMissed, allExtra, allRiseRMSE, allFallRMSE);

                end        
            end

            allR2toPlot{R2counterToPlot} = allR2;

            timeAtMetricGenScenarioEnd = datetime('now');
            disp(['Time after 1 scenario complete: ', char(timeAtMetricGenScenarioEnd)]);
            
            if startsWith(model_params{j}, 'optresult') && strcmpi(model_names{j}, 'Hill') && contains(model_params{j}, 'r01') && ~contains(model_params{j}, 'widthOptOnDF')
                modName = 'Hill-r01';
            elseif startsWith(model_params{j}, 'optresult') && strcmpi(model_names{j}, 'Hill') && contains(model_params{j}, 'r12') && ~contains(model_params{j}, 'widthOptOnDF')
                modName = 'Hill-r12';
            elseif startsWith(model_params{j}, 'network') && strcmpi(model_names{j}, 'NN') && contains(model_params{j}, 'r01')
                modName = 'NN-r01';
            elseif startsWith(model_params{j}, 'network') && strcmpi(model_names{j}, 'NN') && contains(model_params{j}, 'r12')
                modName = 'NN-r12';
            elseif startsWith(model_params{j}, 'Milloptresult') && strcmpi(model_names{j}, 'MillHill') && contains(model_params{j}, 'r01')
                modName = 'MillHill-r01';
            elseif startsWith(model_params{j}, 'Milloptresult') && strcmpi(model_names{j}, 'MillHill') && contains(model_params{j}, 'r12')
                modName = 'MillHill-r12';
            elseif startsWith(model_params{j}, 'standard') && strcmpi(model_names{j}, 'Hill')
                modName = 'std-Hill';
            elseif startsWith(model_params{j}, 'optresult') && strcmpi(model_names{j}, 'Hill') && contains(model_params{j}, 'r01') && contains(model_params{j}, 'widthOptOnDF') 
                modName = 'Hill-r01-wOnDF';
            elseif startsWith(model_params{j}, 'optresult') && strcmpi(model_names{j}, 'Hill') && contains(model_params{j}, 'r12') && contains(model_params{j}, 'widthOptOnDF') 
                modName = 'Hill-r12-wOnDF';
            end

            ModelNameCol{end+1} = modName;
            BirdCol(end+1) = k;
            MuscleCol(end+1) = l;

            avgRMSE(end+1) = mean(allRMSE);
            stdRMSE(end+1) = std(allRMSE);
            avgCorr(end+1) = meanCorrelationCalculation(allCorr);
            medianR2(end+1) = median(allR2);
            minR2(end+1) = min(allR2);
            maxR2(end+1) = max(allR2);

            missedPcnt = 100 * allMissed ./ allMeasured;
            extraPcnt = 100 * allExtra ./ allPredicted;

            meanMissedPcnt(end+1) = sum(allMissed)/sum(allMeasured)*100;
            meanExtraPcnt(end+1) = sum(allExtra)/sum(allPredicted)*100;

            minMissedPcnt(end+1) = min(missedPcnt);
            maxMissedPcnt(end+1) = max(missedPcnt);
            minExtraPcnt(end+1) = min(extraPcnt);
            maxExtraPcnt(end+1) = max(extraPcnt);

            rmseRise(end+1) = mean(allRiseRMSE);
            rmseFall(end+1) = mean(allFallRMSE);
            
        end
    end
end
timeAtMetricGenEnd = datetime('now');
disp(['End of metric generation: ', char(timeAtMetricGenEnd)]);


SummaryTable = table(ModelNameCol', BirdCol', MuscleCol', ...
    avgRMSE', stdRMSE', medianR2', minR2', maxR2', meanMissedPcnt', ...
    meanExtraPcnt', minMissedPcnt', maxMissedPcnt', minExtraPcnt', maxExtraPcnt', ...
    rmseRise', rmseFall', ...
    'VariableNames', {'ModelName', 'Bird', 'Muscle', 'AvgRMSE', 'StdRMSE', ...
    'MedianR2', 'MinR2', 'MaxR2', 'MeanMissedPcnt', 'MeanExtraPcnt', ...
    'MinMissedPcnt', 'MaxMissedPcnt', 'MinExtraPcnt', 'MaxExtraPcnt', ...
    'RMSErise', 'RMSEfall'});

%% plot coefficient of determination 'histograms'
xPos=[1,1.5,2,2.5];
jitterAmount=0.01;
transparency=0.5;
colors=[
0.502, 0.000, 0.125;
0.200, 0.300, 0.800;
0.75, 0.75, 0.75;
0.400, 0.737, 0.718];

modelNamesTable = {'Hill-r01', 'NN-r01', 'Hill-r12', 'NN-r12'};
titlesSubPlots={'Hill-r01', 'NN-r01', 'Hill-r12', 'NN-r12'};
titlesSb={'(A) Models fitted on trial 1', 'Models fitted on trial 12'};
cnt=1;
figure;
for alpha=1:length(modelNamesTable)
    if alpha==1
        subplot(1,2,1);
    elseif alpha==2 || alpha==4
        xPos=[1.2,1.7,2.2,2.7];
    elseif alpha==3
        subplot(1,2,2);
        xPos=[1,1.5,2,2.5];
        cnt=cnt+1;
    end
    hold on;
    for i=1:(length(bird_names)*length(muscle_names))
        idx=strcmp(SummaryTable.ModelName, modelNamesTable{alpha});
        rowsToGetR2From=find(idx);
        r2Val=allR2toPlot{rowsToGetR2From(i)};
        xJit=xPos(i)+(rand(1,length(r2Val))-0.5)*jitterAmount;

        medianV=median(r2Val);
        distFromMedian = abs(r2Val-medianV);
        stdDist=std(r2Val);
        if stdDist==0
            scaledSize=ones(size(r2Val))*120;
        else           
            scaledSize=120-(distFromMedian/stdDist)*80;
            scaledSize=max(min(scaledSize,120),40);
        end
        leg(i)=scatter(xJit, r2Val, scaledSize, 'filled', 'MarkerFaceColor', colors(i,:), 'MarkerFaceAlpha', transparency);

        med=plot([xPos(i)-0.1, xPos(i)+0.1], [medianV, medianV], 'k-', 'LineWidth', 2);
    end
    xlim([0.5 3]);
    ylim([-4 1]);
    if alpha==2 || alpha==4
        xticks([1,1.2,1.5,1.7,2,2.2,2.5,2.7]);
        xticklabels({modelNamesTable{alpha-1},modelNamesTable{alpha},modelNamesTable{alpha-1},modelNamesTable{alpha},modelNamesTable{alpha-1},modelNamesTable{alpha},modelNamesTable{alpha-1},modelNamesTable{alpha}});
    end
    ylabel('R^2');
    title(titlesSb{cnt});
end

l=legend([leg(1), leg(2), leg(3), leg(4), med], {'Bird 1 - LG','Bird 1 - DF','Bird 2 - LG','Bird 2 - DF', 'Median'});
set(l, 'Box', 'off');


%% Force time-series Plots
trial_names = {'Bl3d2_r10_3p8_7cm_Cal', 'Bl3d2_r10_3p8_7cm_Cal', 'BL4d2_r09_3p8_7cm_Cal', 'BL4d2_r09_3p8_7cm_Cal'};
i1 = [r10_ind r10_ind r09_ind r09_ind];
k1 = [1 1 2 2];
l1 = [1 2 1 2];
subPtitles={'(A) Same bird (bird 1), same muscle (LG)', '(B) Same bird (bird 1), different muscle (DF)', '(C) Different bird (bird 2), same muscle (LG)', '(D) Different bird (bird 2), different muscle (DF)'};

figure
colors=lines(8);
for i = 1: length(trial_names)
    subplot(length(trial_names),1,i)
    plot(results(k1(i)).muscles(l1(i)).models(1).trials(i1(i)).time, results(k1(i)).muscles(l1(i)).models(1).trials(i1(i)).measuredForce, 'Color', colors(1,:), 'LineWidth', 1.5, 'DisplayName', 'Measured');
    hold on
    for j = 1:length(modelNamesWTrial)
        if strcmp(modelNamesWTrial{j}, 'NN-r01')
            best=7; %make sure that this index corresponds to the best NN (the one with the lowest objective value) for the NNs that you are using at that moment
            plot(results(k1(i)).muscles(l1(i)).models(j).trials(i1(i)).time, ...
                results(k1(i)).muscles(l1(i)).models(j).trials(i1(i)).nnVariant(best).estimatedForce*results(k1(i)).muscles(l1(i)).models(j).trials(i1(i)).nnVariant(best).fMax, ...
                'Color', colors(j+1,:), 'LineWidth',1.2, 'DisplayName',modelNamesWTrial{j});
        elseif strcmp(modelNamesWTrial{j}, 'NN-r12')
            best=9;
            plot(results(k1(i)).muscles(l1(i)).models(j).trials(i1(i)).time, ...
                results(k1(i)).muscles(l1(i)).models(j).trials(i1(i)).nnVariant(best).estimatedForce*results(k1(i)).muscles(l1(i)).models(j).trials(i1(i)).nnVariant(best).fMax, ...
                'Color', colors(j+1,:), 'LineWidth',1.2, 'DisplayName',modelNamesWTrial{j});
        else            
            plot(results(k1(i)).muscles(l1(i)).models(j).trials(i1(i)).time, ...
                results(k1(i)).muscles(l1(i)).models(j).trials(i1(i)).modelVariant.estimatedForce'*results(k1(i)).muscles(l1(i)).models(j).trials(i1(i)).modelVariant.fMax, ...
                '-.', 'Color', colors(j+1,:), 'LineWidth', 1.2,'DisplayName',modelNamesWTrial{j});
        end
    end
    xlim([7.4 8.9])
    if i==1
        legend('show', 'Location','northeastoutside', 'Box', 'off')
    end
    title(subPtitles{i});
end
xlabel('Time (s)')
ylabel('Force (N)')




function [peakVal, peakLoc, adjDist] = peakDetection(force, muscleName, prescrDist)
    baselineF=median(force);
    center = force-baselineF;
    threshold=mean(center)+1.5*std(center);
    
    [measPk1, measLoc1]=findpeaks(force, 'MinPeakHeight', threshold);
        
    
    if strcmpi(muscleName, 'df') && prescrDist~=0
        adjDist=prescrDist;
    else
        peakDistances=diff(measLoc1);
        avgPeakDistance=mean(peakDistances);
        medPeakDistance=median(peakDistances);
        if medPeakDistance<prescrDist*0.5
            adjDist=prescrDist;
        else
            adjDist=0.7*medPeakDistance;
        end
    end

    
    [peakVal, peakLoc]=findpeaks(force, 'MinPeakHeight', threshold, 'MinPeakDistance', adjDist);
end




function [riseTimes, fallTimes] = computeRiseFallTimes(force, timeArr, peakValues, peakLocations)
    
    meanPeakForce=mean(peakValues);
    halfMeanPeakForce=meanPeakForce*0.5;
    riseTimes=zeros(size(peakLocations));
    fallTimes=zeros(size(peakLocations));

    for ind1=1:length(peakLocations)
        peakIdx=peakLocations(ind1);

        %%Rise time
        forceIdx=peakIdx;
        while forceIdx>1 && force(forceIdx)>halfMeanPeakForce
            forceIdx=forceIdx-1;
        end

        if peakValues(ind1)>halfMeanPeakForce
            if forceIdx>1
                riseTimes(ind1) = timeArr(peakIdx)-timeArr(forceIdx);
            else
                riseTimes(ind1)=NaN;
                
            end
        else
            riseTimes(ind1)=NaN;
            
        end

        %%Fall time
        forceIdx=peakIdx;
        while forceIdx<length(force) && force(forceIdx)>halfMeanPeakForce
            forceIdx=forceIdx+1;
        end

        if peakValues(ind1)>halfMeanPeakForce
            if forceIdx<length(force)
                fallTimes(ind1) = timeArr(forceIdx)-timeArr(peakIdx);
            else
                fallTimes(ind1)=NaN;
                
            end
        else
            fallTimes(ind1)=NaN;
            
        end

    end
end

function [matched, unmatchedMeas, unmatchedPred] = peakMatching(measPeakLoc, predPeakLoc)
    if length(measPeakLoc)<2
        window=0;
        return;
    end

    interv=diff(measPeakLoc);
    robustInterv=prctile(interv, 25);
    window=round(0.45*robustInterv);

    matched=[];
    unmatchedMeas = measPeakLoc;
    unmatchedPred = predPeakLoc;

    usedPreds=false(size(predPeakLoc));

    for i=1:length(measPeakLoc)
        m=measPeakLoc(i);

        idxInWindow=find(abs(predPeakLoc-m)<=window);

        if isempty(idxInWindow)
            continue;
        end

        %if there are more than one peaks within the window, choose the
        %closest predicted one for the measured peak
        validPreds=predPeakLoc(idxInWindow);
        [minDiff, relIdx]=min(abs(validPreds-m));
        predIdx = idxInWindow(relIdx);

        if ~usedPreds(predIdx)
            matched(end+1,:)=[m, predPeakLoc(predIdx), minDiff];
            usedPreds(predIdx)=true;
        end

    end
    matchedMeas = matched(:,1);
    matchedPred = matched(:,2);

    unmatchedMeas = setdiff(measPeakLoc, matchedMeas, 'stable');
    unmatchedPred = setdiff(predPeakLoc, matchedPred, 'stable');
end

function meanCorr=meanCorrelationCalculation(correlationArray)
    z=atanh(correlationArray);
    meanZ=mean(z);
    meanCorr=tanh(meanZ);
end


function [allRMSE, allCorr, allR2, allMeasured, ...
    allPredicted, allMissed, allExtra, allRiseRMSE, ...
    allFallRMSE] = getMetrics(indI, indJ, indK, indL, ...
    force, f_max, innerLoopInd, muscle_names, time, ...
    allRMSE, allCorr, allR2, allMeasured, allPredicted, allMissed, ...
    allExtra, allRiseRMSE, allFallRMSE)

    rmseForceTemp(innerLoopInd)=rmse(force(indI,indJ,indK,indL).est(:), force(indI,1,indK,indL).meas(:)/f_max(indI,indJ,indK,indL));
    allRMSE(end+1)=rmseForceTemp(innerLoopInd);
    
    corrForceTemp(innerLoopInd) = corr(force(indI,indJ,indK,indL).est(:), force(indI,1,indK,indL).meas(:)/f_max(indI,indJ,indK,indL));
    allCorr(end+1)=corrForceTemp(innerLoopInd);
    
    
    [r2(indI,indJ,indK,indL,innerLoopInd), RMSE1(indI,indJ,indK,indL, innerLoopInd)]= rsquare(force(indI,1,indK,indL).meas(:)/f_max(indI,indJ,indK,indL),force(indI,indJ,indK,indL).est(:)); %value to be used for paper                       
    allR2(end+1) = r2(indI,indJ,indK,indL,innerLoopInd);
    
    
    [measPeakValues, measPeakLocations, peakDistanceMeas]=peakDetection(force(indI,1,indK,indL).meas, muscle_names{indL}, 0);                        
    [measRiseTs, measFallTs] = computeRiseFallTimes(force(indI,1,indK,indL).meas, time(indI,indK,indL).t, measPeakValues, measPeakLocations);
    allMeasured(end+1)=length(measPeakValues);                        
    
    
    [predPeakValues, predPeakLocations, peakDistancePred]=peakDetection(force(indI,indJ,indK,indL).est, muscle_names{indL}, peakDistanceMeas);
    [predRiseTs, predFallTs] = computeRiseFallTimes(force(indI,indJ,indK,indL).est, time(indI,indK,indL).t, predPeakValues, predPeakLocations);
    allPredicted(end+1)=length(predPeakValues);
    
    [matchedPairs, unmatchedMeasPeaks, unmatchedPredPeaks] = peakMatching(measPeakLocations, predPeakLocations);
    
    allMissed(end+1)=length(unmatchedMeasPeaks);
    allExtra(end+1)=length(unmatchedPredPeaks);
    
    allRiseRMSE(end+1)=rmse(nanmean(measRiseTs),nanmean(predRiseTs));
    allFallRMSE(end+1)=rmse(nanmean(measFallTs),nanmean(predFallTs));
end