%The following code performs the analysis/testing of a NN model.
%Testing is performed on the reserved trials of the birds that were 
%used for training the NN and also on all trials of the 
%excluded-from-training bird.

clear all
close all
clc

%% Add to path
addpath(genpath('/path/to/code/')) %Change this to the folder where you downloaded this code

%% Settings
bird_names = {'Bl3', 'BL4', 'Or3', 'Ye3', 'Pu1'}; %Birds in order of 1 to 5
muscle_names = {'LG', 'DF'};

%Write here the location of the data
folder_data = '/path/to/NN_Data/';
folder_networksparams = '/path/to/network/and/params/files/'; %Change this to the folder where the networks and model parameter sets are stored

bird_data = readmatrix([folder_data 'MuscleMorphologyData']);
warning('The columns of the xlsx file are hard coded, please make sure that they match your version of MuscleMorphologyData.xlsx in getMuscleParameters.m')


filenamePatterns={'Bl3*network.mat', 'BL4*network.mat','Or3*network.mat','Pu1*network.mat'};
mdlFileArray={};

for i=1:length(filenamePatterns)
    files=dir(fullfile(folder_networksparams,filenamePatterns{i}));
    fullNames=fullfile({files.folder}, {files.name});
    mdlFileArray = [mdlFileArray, fullNames];
end

mdlFileArray=mdlFileArray(:);

%nnFilesArr contains the filenames of the best big NNs from each network
%type (Bl3_all_all, BL4_all_all, etc). The files follow this naming
%convention: '<nameOfExcludedBird>_<obstacleHeight>_<speed>_network'
nnFilesArr = {'Bl3_all_all_4network', 'BL4_all_all_4network', 'Or3_all_all_2network', 'Pu1_all_all_5network'};
modelNamesWBird={'NN-b1', 'NN-b2','NN-b3','NN-b5'};


%% Uncomment to plot NN losses, FLV relationships and training/validation losses of all 'best' NNs
% Plot NN Losses
% plotBigNNLoss();
% 
%Plot FLV relationships
% plotBigNNFLVs();
% 
% % Plot losses of all NN models
% % Comment-out the following for-loop if NN losses are not needed.
% figure;
% tiledlayout(2,2);
% 
% for i=1:length(nnFilesArr)
%     nexttile
%     Mdl = load(nnFilesArr{i});
%     iteration = Mdl.NN.Mdl.TrainingHistory.Iteration;
%     trainLosses = Mdl.NN.Mdl.TrainingHistory.TrainingLoss;
%     valLosses = Mdl.NN.Mdl.TrainingHistory.ValidationLoss;
%     plot(iteration,trainLosses,iteration,valLosses)
%     
%     legend(["Training","Validation"])
%     xlabel("Iteration")
%     ylabel("Mean Squared Error")
% 
%     title(modelNamesWBird{i});
% 
%     %savefig(strcat(nnFilesArr{i},'.fig'));
% end




mdlVariantCounter=1;


ModelNameCol = [];
BirdCol = [];
MuscleCol = [];

%Metric Columns
avgRMSE = [];
stdRMSE = [];
avgRMSEpcnt = [];
stdRMSEpcnt = [];
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


%% Select data for analysis depending on the selected NN model
for k = 1:length(bird_names)
    trial_names = dir([folder_data filesep bird_names{k}]);
    results(k).birdID=bird_names{k};

    %% Constructing the testing dataset depending on which bird was excluded and for each muscle
    for l = 1:length(muscle_names)
        results(k).muscles(l).muscleName=muscle_names{l};
        resultsStructIdx=1;
        for mdlIndx = 1:length(mdlFileArray)
            results(k).muscles(l).models(resultsStructIdx).modelName=...
                modelNamesWBird{resultsStructIdx};
            


            load(mdlFileArray{mdlIndx}); 
            idxToCutChar=find(mdlFileArray{mdlIndx}=='/', 1, 'last');
            newStr=mdlFileArray{mdlIndx}(idxToCutChar+1:end);
            model_params = newStr;
            model_names = 'NN';
            nnFile=mdlFileArray{mdlIndx};

            
            paramStruct = pStruct;
            excludedBird = paramStruct.excludedBird;
            obstacleHeightSetting = paramStruct.obstacleHeight;
            speedSetting = paramStruct.speed;
            testTrialNames = paramStruct.testTrials;

            if contains(mdlFileArray{mdlIndx}, 'Bl3') && mdlIndx==1
                mdlVariantCounter=1;
                R2counterToPlot=R2counterToPlot+1;
                allRMSE = [];
                allR2 = [];
                allCorr = [];
                allMissed = [];
                allExtra = [];
                allMeasured = [];
                allPredicted = [];
                allRiseRMSE = [];
                allFallRMSE = [];
                allRMSEpcnt = [];
            elseif contains(mdlFileArray{mdlIndx}, 'BL4') && mdlIndx==6
                mdlVariantCounter=1;
                R2counterToPlot=R2counterToPlot+1;
                allRMSE = [];
                allR2 = [];
                allCorr = [];
                allMissed = [];
                allExtra = [];
                allMeasured = [];
                allPredicted = [];
                allRiseRMSE = [];
                allFallRMSE = [];
                allRMSEpcnt = [];
            elseif contains(mdlFileArray{mdlIndx}, 'Or3') && mdlIndx==11
                mdlVariantCounter=1;
                R2counterToPlot=R2counterToPlot+1;
                allRMSE = [];
                allR2 = [];
                allCorr = [];
                allMissed = [];
                allExtra = [];
                allMeasured = [];
                allPredicted = [];
                allRiseRMSE = [];
                allFallRMSE = [];
                allRMSEpcnt = [];
            elseif contains(mdlFileArray{mdlIndx}, 'Pu1') && mdlIndx==16
                mdlVariantCounter=1;
                R2counterToPlot=R2counterToPlot+1;
                allRMSE = [];
                allR2 = [];
                allCorr = [];
                allMissed = [];
                allExtra = [];
                allMeasured = [];
                allPredicted = [];
                allRiseRMSE = [];
                allFallRMSE = [];
                allRMSEpcnt = [];
            end

            if strcmp(bird_names{k},excludedBird)
                trial_namesFin = selectTrials(speedSetting, obstacleHeightSetting, trial_names);
            else
                if strcmp(muscle_names{l}, 'LG')
                    for i=1:length(testTrialNames)
                        if contains(testTrialNames{i}, bird_names{k})
                            trial_namesFin = [];
                            trial_namesFin.name = testTrialNames{i};
                        elseif contains(trial_names(i).name, 'Ye3')
                            trial_namesFin=trial_names;                        
                        else
                            continue
                        end
                    end
                else %assume DF muscle
                    trial_namesFin = selectTrials(speedSetting, obstacleHeightSetting, trial_names);
                end
            end
    
            toRemove = [];

            %% Cleaning test trials list from junk or subpar entries
            if length(trial_namesFin) < 3
                
            elseif length(trial_namesFin) >= 3
                clear namesTr
                for iT = 1:length(trial_namesFin)
                    namesTr{iT} = trial_namesFin(iT).name;
                end
                for tempInd = 1:length(trial_namesFin)
                    if strcmp(trial_namesFin(tempInd).name, '.') || strcmp(trial_namesFin(tempInd).name, '..')
                        toRemove = [toRemove, tempInd];
                    end
                    if contains(trial_namesFin(tempInd).name, 'Ye3d1', 'IgnoreCase',true)
                        toRemove = [toRemove, tempInd];
                    end
                end
                trial_namesFin(toRemove) = [];
            end

            for i = 1:length(trial_namesFin)
                %% Removing subpar data from the testing set
                if ~strcmp(trial_namesFin(i).name(end-2:end), 'mat')
                    continue
                end 
                if strcmpi(bird_names{k}, 'ye3') && strcmpi(muscle_names{l}, 'df') %do not use day 1 for DF of ye3
                    if strcmpi(trial_namesFin(i).name(1:5), 'Ye3d1')
                        display(['Removing trial: ', trial_namesFin(i).name]);
                        continue
                    end
                end
                if strcmpi(bird_names{k}, 'ye3') && strcmpi(muscle_names{l}, 'lg') %do not use LG of ye3
                    display(['Removing trial: ', trial_namesFin(i).name, ' Ye3-LG']);
                    continue
                end
                if strcmpi(bird_names{k}, 'pu1') && strcmpi(muscle_names{l}, 'df') %do not use DF of pu1
                    display(['Removing trial: ', trial_namesFin(i).name, ' Pu1-DF']);
                    continue
                end
            
                %% Selecting the index of the test trials whose predictions will be plotted, depending on which NN is being analyzed
                if l == 1
                    clear namesTr
                    for iT = 1:length(trial_namesFin)
                        namesTr{iT} = trial_namesFin(iT).name;
                    end
                    if contains(model_params, 'Pu1_level_all_network')
                        if strcmp(trial_namesFin(i).name, 'Or3d1_r03_1p8_Lev_Cal.mat')
                            Or3_indLG = find(~contains(namesTr, 'r03') == 0);
                        elseif strcmp(trial_namesFin(i).name, 'Pu1d1_r01_1p8_Lev_Cal.mat')
                            Pu1_indLG = find(~contains(namesTr, 'r01') == 0);
                        end
                    elseif contains(model_params, 'Pu1_all_all_5network')
                        if strcmp(trial_namesFin(i).name, 'Or3d1_r05_1p8_7cm_Cal.mat')
                            Or3_indLG = find(~contains(namesTr, 'r05') == 0);
                            
                        elseif strcmp(trial_namesFin(i).name, 'Pu1d1_r09_1p8_7cm_Cal.mat')
                            Pu1_indLG = find(~contains(namesTr, 'r09') == 0);
                            
                        end
                    elseif contains(model_params, 'Or3_level_all_network')
                        if strcmp(trial_namesFin(i).name, 'Or3d1_r01_3p5_Lev_Cal.mat')
                            Or3_indLG = find(~contains(namesTr, 'r01') == 0);
                        elseif strcmp(trial_namesFin(i).name, 'BL4d2_r02_3p8_Lev_Cal.mat')
                            BL4_indLG = find(~contains(namesTr, 'r02') == 0);
                        end
                    elseif contains(model_params, 'Or3_all_all_2network')
                        if strcmp(trial_namesFin(i).name, 'Or3d1_r01_3p5_Lev_Cal.mat')
                            Or3_indLG = find(~contains(namesTr, 'r01') == 0);
                            
                        elseif strcmp(trial_namesFin(i).name, 'BL4d2_r02_3p8_Lev_Cal.mat')
                            BL4_indLG = find(~contains(namesTr, 'r02') == 0);
                            
                        end
                    elseif contains(model_params, 'BL4_level_all_network')
                        if strcmp(trial_namesFin(i).name, 'BL4d2_r01_1p8_Lev_Cal.mat')
                            BL4_indLG = find(~contains(namesTr, 'r01') == 0);
                        elseif strcmp(trial_namesFin(i).name, 'Or3d1_r03_1p8_Lev_Cal.mat')
                            Or3_indLG = find(~contains(namesTr, 'r03') == 0);
                        end
                    elseif contains(model_params, 'BL4_all_all_4network')
                        if strcmp(trial_namesFin(i).name, 'Or3d1_r05_1p8_7cm_Cal.mat')
                            Or3_indLG = find(~contains(namesTr, 'r05') == 0);
                            
                        elseif strcmp(trial_namesFin(i).name, 'BL4d2_r05_1p8_7cm_Cal.mat')
                            BL4_indLG = find(~contains(namesTr, 'r05') == 0);
                            
                        end
                    elseif contains(model_params, 'Bl3_level_all_network')
                        if strcmp(trial_namesFin(i).name, 'Bl3d2_r01_1p8_Lev_Cal.mat')
                            Bl3_indLG = find(~contains(namesTr, 'r01') == 0);
                        elseif strcmp(trial_namesFin(i).name, 'Or3d1_r03_1p8_Lev_Cal.mat')
                            Or3_indLG = find(~contains(namesTr, 'r03') == 0);
                        end
                    elseif contains(model_params, 'Bl3_all_all_4network')
                        if strcmp(trial_namesFin(i).name, 'Or3d1_r05_1p8_7cm_Cal.mat')
                            Or3_indLG = find(~contains(namesTr, 'r05') == 0);
                            
                        elseif strcmp(trial_namesFin(i).name, 'Bl3d2_r11_1p8_7cm_Cal.mat')
                            Bl3_indLG = find(~contains(namesTr, 'r11') == 0);
                            
                        end
                    end
                        
                elseif l==2
                    clear namesTr
                    for iT = 1:length(trial_namesFin)
                        namesTr{iT} = trial_namesFin(iT).name;
                    end
                    if contains(model_params, 'Pu1_level_all_network')
                        if strcmp(trial_namesFin(i).name, 'Or3d1_r03_1p8_Lev_Cal.mat')
                            Or3_indDF = find(~contains(namesTr, 'r03') == 0);
                        elseif strcmp(trial_namesFin(i).name, 'Ye3d2_t01_1p8_Lev_Cal.mat')
                            Ye3_indDF = find(~contains(namesTr, 'd2_t01') == 0);
                        end
                    elseif contains(model_params, 'Pu1_all_all_5network')
                        if strcmp(trial_namesFin(i).name, 'Or3d1_r05_1p8_7cm_Cal.mat')
                            Or3_indDF = find(~contains(namesTr, 'r05') == 0);
                            
                        elseif strcmp(trial_namesFin(i).name, 'Ye3d2_t12_1p8_7cm_Cal.mat')
                            Ye3_indDF = find(~contains(namesTr, 'd2_t12') == 0);
                            
                        end
                    elseif contains(model_params, 'Or3_level_all_network')
                        if strcmp(trial_namesFin(i).name, 'Or3d1_r01_3p5_Lev_Cal.mat')
                            Or3_indDF = find(~contains(namesTr, 'r01') == 0);
                        elseif strcmp(trial_namesFin(i).name, 'BL4d2_r02_3p8_Lev_Cal.mat')
                            BL4_indDF = find(~contains(namesTr, 'r02') == 0);
                        end
                    elseif contains(model_params, 'Or3_all_all_2network')
                        if strcmp(trial_namesFin(i).name, 'Or3d1_r01_3p5_Lev_Cal.mat')
                            Or3_indDF = find(~contains(namesTr, 'r01') == 0);
                            
                        elseif strcmp(trial_namesFin(i).name, 'BL4d2_r02_3p8_Lev_Cal.mat')
                            BL4_indDF = find(~contains(namesTr, 'r02') == 0);
                            
                        end
                    elseif contains(model_params, 'BL4_level_all_network')
                        if strcmp(trial_namesFin(i).name, 'BL4d2_r01_1p8_Lev_Cal.mat')
                            BL4_indDF = find(~contains(namesTr, 'r01') == 0);
                        elseif strcmp(trial_namesFin(i).name, 'Or3d1_r03_1p8_Lev_Cal.mat')
                            Or3_indDF = find(~contains(namesTr, 'r03') == 0);
                        end
                    elseif contains(model_params, 'BL4_all_all_4network')
                        if strcmp(trial_namesFin(i).name, 'Or3d1_r05_1p8_7cm_Cal.mat')
                            Or3_indDF = find(~contains(namesTr, 'r05') == 0);
                            
                        elseif strcmp(trial_namesFin(i).name, 'BL4d2_r05_1p8_7cm_Cal.mat')
                            BL4_indDF = find(~contains(namesTr, 'r05') == 0);
                            
                        end
                    elseif contains(model_params, 'Bl3_level_all_network')
                        if strcmp(trial_namesFin(i).name, 'Bl3d2_r01_1p8_Lev_Cal.mat')
                            Bl3_indDF = find(~contains(namesTr, 'r01') == 0);
                        elseif strcmp(trial_namesFin(i).name, 'Or3d1_r03_1p8_Lev_Cal.mat')
                            Or3_indDF = find(~contains(namesTr, 'r03') == 0);
                        end
                    elseif contains(model_params, 'Bl3_all_all_4network')
                        if strcmp(trial_namesFin(i).name, 'Or3d1_r05_1p8_7cm_Cal.mat')
                            Or3_indDF = find(~contains(namesTr, 'r05') == 0);
                            
                        elseif strcmp(trial_namesFin(i).name, 'Bl3d2_r11_1p8_7cm_Cal.mat')
                            Bl3_indDF = find(~contains(namesTr, 'r11') == 0);
                            
                        end
                    end
                end
    
                % Indices explanation: i=trial name, k=bird name/number, 
                % l=muscle name/number, j=network ID
                %% Get forces, calculate RMSEs and correlations
                [time(i,k,l).t, data(i,k,l).l_ce, data(i,k,l).v_ce, data(i,k,l).EMG, force(i,1,k,l).meas] = loadDataFile(bird_names{k}, muscle_names{l}, trial_namesFin(i).name,folder_data); 
                results(k).muscles(l).models(resultsStructIdx).trials(i).measuredForce=force(i,1,k,l).meas;
                results(k).muscles(l).models(resultsStructIdx).trials(i).time=time(i,k,l).t;
                
               
                [force(i,mdlIndx,k,l).est,f_max(i,mdlIndx,k,l)] = getForce(folder_networksparams, data(i,k,l).l_ce, ...
                                                                  data(i,k,l).v_ce, data(i,k,l).EMG, model_names, ...
                                                                  model_params, muscle_names{l}, bird_names{k}, ...
                                                                  bird_data);

                results(k).muscles(l).models(resultsStructIdx).trials(i).nnVariant(mdlVariantCounter).nnVariantID=mdlVariantCounter;
                results(k).muscles(l).models(resultsStructIdx).trials(i).nnVariant(mdlVariantCounter).estimatedForce=force(i,mdlIndx,k,l).est;
                results(k).muscles(l).models(resultsStructIdx).trials(i).nnVariant(mdlVariantCounter).fMax=f_max(i,mdlIndx,k,l);

                [allRMSE, allCorr, allR2, allMeasured, allPredicted, allMissed, allExtra, allRiseRMSE, allFallRMSE, allRMSEpcnt] = ...
                                                                    getMetrics(i, mdlIndx, k, l, force, f_max, i, muscle_names, ...
                                                                    time, allRMSE, allCorr, allR2, allMeasured, allPredicted, ...
                                                                    allMissed, allExtra, allRiseRMSE, allFallRMSE, allRMSEpcnt);
            end
            mdlVariantCounter=mdlVariantCounter+1;

            if strcmpi(bird_names{k}, 'Ye3') && strcmpi(muscle_names{l}, 'LG')
                continue
            elseif strcmpi(bird_names{k}, 'Pu1') && strcmpi(muscle_names{l}, 'DF')
                continue
            end

            if contains(mdlFileArray{mdlIndx}, 'Bl3') && mdlIndx==5
                resultsStructIdx=resultsStructIdx+1;
                modName = 'NN-b1';
                allR2toPlot{R2counterToPlot} = allR2;
                [ModelNameCol, BirdCol, MuscleCol, avgRMSE, stdRMSE, avgRMSEpcnt, stdRMSEpcnt, avgCorr,...
                    medianR2, minR2, maxR2, meanMissedPcnt, meanExtraPcnt, minMissedPcnt,...
                    maxMissedPcnt, minExtraPcnt, maxExtraPcnt, rmseRise, rmseFall] = ...
                                    perModelTypeMetric(ModelNameCol, BirdCol, MuscleCol, avgRMSE, ...
                                    stdRMSE, avgRMSEpcnt, stdRMSEpcnt, avgCorr, medianR2, minR2, maxR2, meanMissedPcnt, ...
                                    meanExtraPcnt, minMissedPcnt,maxMissedPcnt, minExtraPcnt, ...
                                    maxExtraPcnt, rmseRise, rmseFall,...
                                    modName, k, l, allRMSE, allCorr, allR2, ...
                                    allMissed, allMeasured, allExtra, allPredicted, ...
                                    allRiseRMSE, allFallRMSE, allRMSEpcnt);                

            elseif contains(mdlFileArray{mdlIndx}, 'BL4') && mdlIndx==10
                resultsStructIdx=resultsStructIdx+1;
                modName = 'NN-b2';
                allR2toPlot{R2counterToPlot} = allR2;
                [ModelNameCol, BirdCol, MuscleCol, avgRMSE, stdRMSE, avgRMSEpcnt, stdRMSEpcnt, avgCorr,...
                    medianR2, minR2, maxR2, meanMissedPcnt, meanExtraPcnt, minMissedPcnt,...
                    maxMissedPcnt, minExtraPcnt, maxExtraPcnt, rmseRise, rmseFall] = ...
                                    perModelTypeMetric(ModelNameCol, BirdCol, MuscleCol, avgRMSE, ...
                                    stdRMSE, avgRMSEpcnt, stdRMSEpcnt, avgCorr, medianR2, minR2, maxR2, meanMissedPcnt, ...
                                    meanExtraPcnt, minMissedPcnt,maxMissedPcnt, minExtraPcnt, ...
                                    maxExtraPcnt, rmseRise, rmseFall,...
                                    modName, k, l, allRMSE, allCorr, allR2, ...
                                    allMissed, allMeasured, allExtra, allPredicted, ...
                                    allRiseRMSE, allFallRMSE, allRMSEpcnt);
            elseif contains(mdlFileArray{mdlIndx}, 'Or3') && mdlIndx==15
                resultsStructIdx=resultsStructIdx+1;
                modName = 'NN-b3';
                allR2toPlot{R2counterToPlot} = allR2;
                [ModelNameCol, BirdCol, MuscleCol, avgRMSE, stdRMSE, avgRMSEpcnt, stdRMSEpcnt, avgCorr,...
                    medianR2, minR2, maxR2, meanMissedPcnt, meanExtraPcnt, minMissedPcnt,...
                    maxMissedPcnt, minExtraPcnt, maxExtraPcnt, rmseRise, rmseFall] = ...
                                    perModelTypeMetric(ModelNameCol, BirdCol, MuscleCol, avgRMSE, ...
                                    stdRMSE, avgRMSEpcnt, stdRMSEpcnt, avgCorr, medianR2, minR2, maxR2, meanMissedPcnt, ...
                                    meanExtraPcnt, minMissedPcnt,maxMissedPcnt, minExtraPcnt, ...
                                    maxExtraPcnt, rmseRise, rmseFall,...
                                    modName, k, l, allRMSE, allCorr, allR2, ...
                                    allMissed, allMeasured, allExtra, allPredicted, ...
                                    allRiseRMSE, allFallRMSE, allRMSEpcnt);
                
            elseif contains(mdlFileArray{mdlIndx}, 'Pu1') && mdlIndx==20
                resultsStructIdx=resultsStructIdx+1;
                modName = 'NN-b5';
                allR2toPlot{R2counterToPlot} = allR2;
                [ModelNameCol, BirdCol, MuscleCol, avgRMSE, stdRMSE, avgRMSEpcnt, stdRMSEpcnt, avgCorr,...
                    medianR2, minR2, maxR2, meanMissedPcnt, meanExtraPcnt, minMissedPcnt,...
                    maxMissedPcnt, minExtraPcnt, maxExtraPcnt, rmseRise, rmseFall] = ...
                                    perModelTypeMetric(ModelNameCol, BirdCol, MuscleCol, avgRMSE, ...
                                    stdRMSE, avgRMSEpcnt, stdRMSEpcnt, avgCorr, medianR2, minR2, maxR2, meanMissedPcnt, ...
                                    meanExtraPcnt, minMissedPcnt,maxMissedPcnt, minExtraPcnt, ...
                                    maxExtraPcnt, rmseRise, rmseFall,...
                                    modName, k, l, allRMSE, allCorr, allR2, ...
                                    allMissed, allMeasured, allExtra, allPredicted, ...
                                    allRiseRMSE, allFallRMSE, allRMSEpcnt);
            end
        end
    end
end


SummaryTable = table(ModelNameCol', BirdCol', MuscleCol', ...
    avgRMSE', stdRMSE', avgRMSEpcnt', stdRMSEpcnt', medianR2', minR2', maxR2', meanMissedPcnt', ...
    meanExtraPcnt', minMissedPcnt', maxMissedPcnt', minExtraPcnt', maxExtraPcnt', ...
    rmseRise', rmseFall', ...
    'VariableNames', {'ModelName', 'Bird', 'Muscle', 'AvgRMSE', 'StdRMSE', 'AvgRMSEpcnt', 'StdRMSEpcnt', ...
    'MedianR2', 'MinR2', 'MaxR2', 'MeanMissedPcnt', 'MeanExtraPcnt', ...
    'MinMissedPcnt', 'MaxMissedPcnt', 'MinExtraPcnt', 'MaxExtraPcnt', ...
    'RMSErise', 'RMSEfall'});



%% plot coefficient of determination 'histograms'
allR2toPlot2=allR2toPlot(~cellfun('isempty',allR2toPlot));

xPos1=[1,2,3,4,5,6,7,8,9,10];
xPos2=[1.2,2.2,3.2,4.2,5.2,6.2,7.2,8.2,9.2,10.2];
xPos3=[1.4,2.4,3.4,4.4,5.4,6.4,7.4,8.4,9.4,10.4];
xPos4=[1.6,2.6,3.6,4.6,5.6,6.6,7.6,8.6,9.6,10.6];
jitterAmount=0.01;
transparency=0.5;
colors=[
0.502, 0.000, 0.125;
0.200, 0.300, 0.800;
0.75, 0.75, 0.75;
0.400, 0.737, 0.718;
0.850,0.325,0.098;
0.929,0.694,0.125;
0.494,0.184,0.556;
0.466,0.674,0.188;
0.301,0.745,0.933;
0.172,0.627,0.172];

modelNamesTable = {'NN-b1','NN-b2','NN-b3','NN-b5'};
titlesSubPlots={'NN-b1','NN-b2','NN-b3','NN-b5'};
figure;
for alpha=1:length(modelNamesTable)
    if alpha==1
        xPos=xPos1;
    elseif alpha==2 
        xPos=xPos2;
    elseif alpha==3       
        xPos=xPos3;
    elseif alpha==4       
        xPos=xPos4;
    end
    hold on;
    for i=1:8
        idx=strcmp(SummaryTable.ModelName, modelNamesTable{alpha});
        rowsToGetR2From=find(idx);
        r2Val=allR2toPlot2{rowsToGetR2From(i)};
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

    if alpha==1
        xticks([1,1.2,1.4,1.6,2,2.2,2.4,2.6,3,3.2,3.4,3.6,4,4.2,4.4,4.6,5,5.2,5.4,5.6,6,6.2,6.4,6.6,7,7.2,7.4,7.6,8,8.2,8.4,8.6]);
        xticklabels({modelNamesTable{alpha},modelNamesTable{alpha+1},modelNamesTable{alpha+2},modelNamesTable{alpha+3}, ...
                    modelNamesTable{alpha},modelNamesTable{alpha+1},modelNamesTable{alpha+2},modelNamesTable{alpha+3}, ...
                    modelNamesTable{alpha},modelNamesTable{alpha+1},modelNamesTable{alpha+2},modelNamesTable{alpha+3}, ...
                    modelNamesTable{alpha},modelNamesTable{alpha+1},modelNamesTable{alpha+2},modelNamesTable{alpha+3}, ...
                    modelNamesTable{alpha},modelNamesTable{alpha+1},modelNamesTable{alpha+2},modelNamesTable{alpha+3}, ...
                    modelNamesTable{alpha},modelNamesTable{alpha+1},modelNamesTable{alpha+2},modelNamesTable{alpha+3}, ...
                    modelNamesTable{alpha},modelNamesTable{alpha+1},modelNamesTable{alpha+2},modelNamesTable{alpha+3}, ...
                    modelNamesTable{alpha},modelNamesTable{alpha+1},modelNamesTable{alpha+2},modelNamesTable{alpha+3}});
    end
    ylabel('R^2');
end
title('R^2 Distribution of all Models across all Test Scenarios');
l=legend([leg(1), leg(2), leg(3), leg(4),leg(5),leg(6),leg(7),leg(8), med], ...
    {'Bird 1 - LG','Bird 1 - DF','Bird 2 - LG','Bird 2 - DF', ...
    'Bird 3 - LG','Bird 3 - DF','Bird 4 - DF','Bird 5 - LG', ...
    'Median'});
set(l, 'Box', 'off');




%% Force Plots
for alpha=1:length(nnFilesArr)
    
    if strcmp(nnFilesArr{alpha}, 'Pu1_level_all_network') || strcmp(nnFilesArr{alpha}, 'Pu1_all_all_5network')
        varIdx=5;
        mdlIdx=4;
        i1 = [1 5 9 12];
        k1 = [3 3 5 4];    
    elseif strcmp(nnFilesArr{alpha}, 'Or3_level_all_network') || strcmp(nnFilesArr{alpha}, 'Or3_all_all_2network')
        varIdx=2;
        mdlIdx=3;
        i1 = [1 1 1 2];
        k1 = [3 3 2 2];
    elseif strcmp(nnFilesArr{alpha}, 'BL4_level_all_network') || strcmp(nnFilesArr{alpha}, 'BL4_all_all_4network')
        varIdx=4;
        mdlIdx=2;
        i1 = [1 5 5 5];
        k1 = [3 3 2 2];
    elseif strcmp(nnFilesArr{alpha}, 'Bl3_level_all_network') || strcmp(nnFilesArr{alpha}, 'Bl3_all_all_4network')
        varIdx=4;
        mdlIdx=1;
        i1 = [1 5 11 11];
        k1 = [3 3 1 1];
    end
    
    l1 = [1 2 1 2]; 
    
    figure
    for i = 1: length(i1)
        subplot(length(i1),1,i)
    
        plot(results(k1(i)).muscles(l1(i)).models(1).trials(i1(i)).time, results(k1(i)).muscles(l1(i)).models(1).trials(i1(i)).measuredForce, 'LineWidth', 2, 'DisplayName', 'Measured');
        hold on
        for j = 1:length(model_names)
            plot(results(k1(i)).muscles(l1(i)).models(mdlIdx).trials(i1(i)).time, results(k1(i)).muscles(l1(i)).models(mdlIdx).trials(i1(i)).nnVariant(varIdx).estimatedForce, 'LineWidth', 2);
              
        end
        xlim([10 12])
    end
    legend({'Measurements', model_names{:}})
end



























function trials = selectTrials(speedSetting, obstacleHeightSetting, trialsIn)
    if strcmpi(speedSetting, 'all') && strcmpi(obstacleHeightSetting, 'all')
        trials = trialsIn;
    elseif strcmpi(speedSetting, '1p8') && strcmpi(obstacleHeightSetting, 'all')
        substrings = {'3p8', '4p5', '3p0', '3p5'};
        mask = contains(trialsIn.name, substrings); %make array of indices for trials including the above substrings in their name
        trialsIn(mask) = [];
        trials=trialsIn;
    elseif strcmpi(speedSetting, 'all') && strcmpi(obstacleHeightSetting, 'level')
        substrings = {'5cm', '7cm'};
        trialsTemp = strings(15,1);
        for i = 1:length(trialsIn)
            trialsTemp(i,1) = trialsIn(i).name;
        end
        mask = contains(trialsTemp, substrings); %make array of indices for trials including the above substrings in their name
        trialsIn(mask) = [];
        trials = trialsIn;
    elseif strcmpi(speedSetting, '1p8') && strcmpi(obstacleHeightSetting, 'level')
        warning("Currently not implemented. Too few training trials remain under this condition, not allowing NN training and testing in a typical ML scenario.");
    end
end



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
    allFallRMSE, allRMSEpcnt] = getMetrics(indI, indJ, indK, indL, ...
    force, f_max, innerLoopInd, muscle_names, time, ...
    allRMSE, allCorr, allR2, allMeasured, allPredicted, allMissed, ...
    allExtra, allRiseRMSE, allFallRMSE, allRMSEpcnt)

    rmseForceTemp(innerLoopInd)=rmse(force(indI,indJ,indK,indL).est(:), force(indI,1,indK,indL).meas(:)/f_max(indI,indJ,indK,indL));
    allRMSE(end+1)=rmseForceTemp(innerLoopInd);
    
    
    maxTrialForceTemp(innerLoopInd) = max(force(indI,1,indK,indL).meas(:)/f_max(indI,indJ,indK,indL));
    RMSEpcntTemp(innerLoopInd) = rmseForceTemp(innerLoopInd)/maxTrialForceTemp(innerLoopInd)*100;
    allRMSEpcnt(end+1)=RMSEpcntTemp(innerLoopInd);


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

function [ModelNameCol, BirdCol, MuscleCol, avgRMSE, stdRMSE, avgRMSEpcnt, stdRMSEpcnt, avgCorr,...
        medianR2, minR2, maxR2, meanMissedPcnt, meanExtraPcnt, minMissedPcnt,...
        maxMissedPcnt, minExtraPcnt, maxExtraPcnt, rmseRise, rmseFall] = ...
                    perModelTypeMetric(ModelNameCol, BirdCol, MuscleCol, avgRMSE, ...
                    stdRMSE, avgRMSEpcnt, stdRMSEpcnt, avgCorr, medianR2, minR2, maxR2, meanMissedPcnt, ...
                    meanExtraPcnt, minMissedPcnt,maxMissedPcnt, minExtraPcnt, ...
                    maxExtraPcnt, rmseRise, rmseFall,...
                    modName, k, l, allRMSE, allCorr, allR2, ...
                    allMissed, allMeasured, allExtra, allPredicted, ...
                    allRiseRMSE, allFallRMSE, allRMSEpcnt)

    ModelNameCol{end+1} = modName;
    BirdCol(end+1) = k;
    MuscleCol(end+1) = l;

    avgRMSE(end+1) = mean(allRMSE);
    stdRMSE(end+1) = std(allRMSE);
    avgRMSEpcnt(end+1) = mean(allRMSEpcnt);
    stdRMSEpcnt(end+1)=std(allRMSEpcnt);
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