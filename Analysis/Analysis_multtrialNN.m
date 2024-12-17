%The following code performs the analysis/testing of a NN model.
%Testing is performed on the reserved trials of the birds that were 
%used for training the NN and also on all trials of the 
%excluded-from-training bird.

clear all
close all
clc

%% Add to path
addpath(genpath('path/to/code/')) %Change this to the folder where you downloaded this code

%% Settings
bird_names = {'Bl3', 'BL4', 'Or3', 'Ye3', 'Pu1'}; %Birds in order of 1 to 5
muscle_names = {'LG', 'DF'};

%Write here the location of the data
folder_base = 'path/to/data/';

bird_data = xlsread([folder_base 'MuscleMorphologyData']);

% Specify in 'nnfile' the network name for which you would like to run the analysis
nnFile = '<nameOfExcludedBird>_<obstacleHeight>_<speed>_network';
nnFilesArr = {'Bl3_all_all_network', 'Pu1_all_all_network', 'Or3_all_all_network', 'BL4_all_all_network'};
load(nnFile); 
model_params = {nnFile};
model_names = {'NN'};

paramStruct = pStruct;
excludedBird = paramStruct.excludedBird;
obstacleHeightSetting = paramStruct.obstacleHeight;
speedSetting = paramStruct.speed;
testTrialNames = paramStruct.testTrials;


%% Plot losses of all NN models
% Comment-out the following for-loop if NN losses are not needed. 
for i=1:length(nnFilesArr)
    Mdl = load(nnFilesArr{i});
    iteration = Mdl.NN.Mdl.TrainingHistory.Iteration;
    trainLosses = Mdl.NN.Mdl.TrainingHistory.TrainingLoss;
    valLosses = Mdl.NN.Mdl.TrainingHistory.ValidationLoss;
    plot(iteration,trainLosses,iteration,valLosses)
    
    legend(["Training","Validation"])
    xlabel("Iteration")
    ylabel("Mean Squared Error")

    savefig(strcat(nnFilesArr{i},'.fig'));
end

%% Select data for analysis depending on the selected NN model
for k = 1:length(bird_names)
    trial_names = dir([folder_base filesep bird_names{k}]);

    %% Constructing the testing dataset depending on which bird was excluded and for each muscle
    for l = 1:length(muscle_names)
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
                display(['Removing trial: ', trial_namesFin(i).name]);
                continue
            end
            if strcmpi(bird_names{k}, 'pu1') && strcmpi(muscle_names{l}, 'df') %do not use DF of pu1
                display(['Removing trial: ', trial_namesFin(i).name]);
                continue
            end
        
            %% Selecting the index of the test trials whose predictions will be plotted, depending on which NN is being analyzed
            if l == 1
                clear namesTr
                for iT = 1:length(trial_namesFin)
                    namesTr{iT} = trial_namesFin(iT).name;
                end
                if strcmp(nnFile, 'Pu1_level_all_network')
                    if strcmp(trial_namesFin(i).name, 'Or3d1_r03_1p8_Lev_Cal.mat')
                        Or3_indLG = find(~contains(namesTr, 'r03') == 0);
                    elseif strcmp(trial_namesFin(i).name, 'Pu1d1_r01_1p8_Lev_Cal.mat')
                        Pu1_indLG = find(~contains(namesTr, 'r01') == 0);
                    end
                elseif strcmp(nnFile, 'Pu1_all_all_network')
                    if strcmp(trial_namesFin(i).name, 'Or3d1_r05_1p8_7cm_Cal.mat')
                        Or3_indLG = find(~contains(namesTr, 'r05') == 0);
                    elseif strcmp(trial_namesFin(i).name, 'Pu1d1_r09_1p8_7cm_Cal.mat')
                        Pu1_indLG = find(~contains(namesTr, 'r09') == 0);
                    end
                elseif strcmp(nnFile, 'Or3_level_all_network')
                    if strcmp(trial_namesFin(i).name, 'Or3d1_r01_3p5_Lev_Cal.mat')
                        Or3_indLG = find(~contains(namesTr, 'r01') == 0);
                    elseif strcmp(trial_namesFin(i).name, 'BL4d2_r02_3p8_Lev_Cal.mat')
                        BL4_indLG = find(~contains(namesTr, 'r02') == 0);
                    end
                elseif strcmp(nnFile, 'Or3_all_all_network')
                    if strcmp(trial_namesFin(i).name, 'Or3d1_r01_3p5_Lev_Cal.mat')
                        Or3_indLG = find(~contains(namesTr, 'r01') == 0);
                    elseif strcmp(trial_namesFin(i).name, 'BL4d2_r02_3p8_Lev_Cal.mat')
                        BL4_indLG = find(~contains(namesTr, 'r02') == 0);
                    end
                elseif strcmp(nnFile, 'BL4_level_all_network')
                    if strcmp(trial_namesFin(i).name, 'BL4d2_r01_1p8_Lev_Cal.mat')
                        BL4_indLG = find(~contains(namesTr, 'r01') == 0);
                    elseif strcmp(trial_namesFin(i).name, 'Or3d1_r03_1p8_Lev_Cal.mat')
                        Or3_indLG = find(~contains(namesTr, 'r03') == 0);
                    end
                elseif strcmp(nnFile, 'BL4_all_all_network')
                    if strcmp(trial_namesFin(i).name, 'Or3d1_r05_1p8_7cm_Cal.mat')
                        Or3_indLG = find(~contains(namesTr, 'r05') == 0);
                    elseif strcmp(trial_namesFin(i).name, 'BL4d2_r05_1p8_7cm_Cal.mat')
                        BL4_indLG = find(~contains(namesTr, 'r05') == 0);
                    end
                elseif strcmp(nnFile, 'Bl3_level_all_network')
                    if strcmp(trial_namesFin(i).name, 'Bl3d2_r01_1p8_Lev_Cal.mat')
                        Bl3_indLG = find(~contains(namesTr, 'r01') == 0);
                    elseif strcmp(trial_namesFin(i).name, 'Or3d1_r03_1p8_Lev_Cal.mat')
                        Or3_indLG = find(~contains(namesTr, 'r03') == 0);
                    end
                elseif strcmp(nnFile, 'Bl3_all_all_network')
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
                if strcmp(nnFile, 'Pu1_level_all_network')
                    if strcmp(trial_namesFin(i).name, 'Or3d1_r03_1p8_Lev_Cal.mat')
                        Or3_indDF = find(~contains(namesTr, 'r03') == 0);
                    elseif strcmp(trial_namesFin(i).name, 'Ye3d2_t01_1p8_Lev_Cal.mat')
                        Ye3_indDF = find(~contains(namesTr, 'd2_t01') == 0);
                    end
                elseif strcmp(nnFile, 'Pu1_all_all_network')
                    if strcmp(trial_namesFin(i).name, 'Or3d1_r05_1p8_7cm_Cal.mat')
                        Or3_indDF = find(~contains(namesTr, 'r05') == 0);
                    elseif strcmp(trial_namesFin(i).name, 'Ye3d2_t12_1p8_7cm_Cal.mat')
                        Ye3_indDF = find(~contains(namesTr, 'd2_t12') == 0);
                    end
                elseif strcmp(nnFile, 'Or3_level_all_network')
                    if strcmp(trial_namesFin(i).name, 'Or3d1_r01_3p5_Lev_Cal.mat')
                        Or3_indDF = find(~contains(namesTr, 'r01') == 0);
                    elseif strcmp(trial_namesFin(i).name, 'BL4d2_r02_3p8_Lev_Cal.mat')
                        BL4_indDF = find(~contains(namesTr, 'r02') == 0);
                    end
                elseif strcmp(nnFile, 'Or3_all_all_network')
                    if strcmp(trial_namesFin(i).name, 'Or3d1_r01_3p5_Lev_Cal.mat')
                        Or3_indDF = find(~contains(namesTr, 'r01') == 0);
                    elseif strcmp(trial_namesFin(i).name, 'BL4d2_r02_3p8_Lev_Cal.mat')
                        BL4_indDF = find(~contains(namesTr, 'r02') == 0);
                    end
                elseif strcmp(nnFile, 'BL4_level_all_network')
                    if strcmp(trial_namesFin(i).name, 'BL4d2_r01_1p8_Lev_Cal.mat')
                        BL4_indDF = find(~contains(namesTr, 'r01') == 0);
                    elseif strcmp(trial_namesFin(i).name, 'Or3d1_r03_1p8_Lev_Cal.mat')
                        Or3_indDF = find(~contains(namesTr, 'r03') == 0);
                    end
                elseif strcmp(nnFile, 'BL4_all_all_network')
                    if strcmp(trial_namesFin(i).name, 'Or3d1_r05_1p8_7cm_Cal.mat')
                        Or3_indDF = find(~contains(namesTr, 'r05') == 0);
                    elseif strcmp(trial_namesFin(i).name, 'BL4d2_r05_1p8_7cm_Cal.mat')
                        BL4_indDF = find(~contains(namesTr, 'r05') == 0);
                    end
                elseif strcmp(nnFile, 'Bl3_level_all_network')
                    if strcmp(trial_namesFin(i).name, 'Bl3d2_r01_1p8_Lev_Cal.mat')
                        Bl3_indDF = find(~contains(namesTr, 'r01') == 0);
                    elseif strcmp(trial_namesFin(i).name, 'Or3d1_r03_1p8_Lev_Cal.mat')
                        Or3_indDF = find(~contains(namesTr, 'r03') == 0);
                    end
                elseif strcmp(nnFile, 'Bl3_all_all_network')
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
            [time(i,k,l).t, data(i,k,l).l_ce, data(i,k,l).v_ce, data(i,k,l).EMG, force(i,1,k,l).meas] = loadDataFile(bird_names{k}, muscle_names{l}, trial_namesFin(i).name,folder_base); 
            for j = 1:length(model_names)
                display(['Predicting force for muscle ', muscle_names{l}, ' and trial: ', trial_namesFin(i).name]);
                [force(i,j,k,l).est,f_max(i,j,k,l)] = getForce(data(i,k,l).l_ce, data(i,k,l).v_ce, data(i,k,l).EMG, model_names{j}, model_params{j}, muscle_names{l}, bird_names{k}, bird_data);
     
                RMSE(i,j,k,l) = rmse(force(i,j,k,l).est(:), force(i,1,k,l).meas(:)/f_max(i,j,k,l));
                maxTrialForce(i,j,k,l) = max(force(i,1,k,l).meas(:)/f_max(i,j,k,l));
                RMSEpcnt(i,j,k,l) = RMSE(i,j,k,l)/maxTrialForce(i,j,k,l)*100; 

                corrP(i,j,k,l) = corr(force(i,j,k,l).est(:), force(i,1,k,l).meas(:)/f_max(i,j,k,l));
                [r2(i,j,k,l), RMSE1(i,j,k,l)] = rsquare(force(i,1,k,l).meas(:)/f_max(i,j,k,l),force(i,j,k,l).est(:));
            end
        end
    end
end

%% Calculate mean RMSE and correlations
for j = 1:length(model_names)
    for k = 1:length(bird_names)
        for l = 1:length(muscle_names)
            RMSE_now = RMSE(:,j,k,l);
            corrP_now = corrP(:, j,k,l);
            RMSE1_now = RMSE1(:,j,k,l);
            R2_now = r2(:,j,k,l);
            RMSEpcnt_now = RMSEpcnt(:,j,k,l);
            

            % identify and remove zeros
            ind_rem = find(abs(RMSE_now)<1e-6);

            RMSE_now(ind_rem) = [];
            RMSE_avg(j,k,l) = mean(RMSE_now);

            corrP_now(ind_rem) = [];
            corrP_avg(j,k,l) = mean(corrP_now);

            RMSE1_now(ind_rem) = [];
            RMSE1_avg(j,k,l) = mean(RMSE1_now);

            R2_now(ind_rem) = [];
            R2_avg(j,k,l) = mean(R2_now);

            RMSEpcnt_now(ind_rem) = [];
            RMSEpcnt_avg(j,k,l) = mean(RMSEpcnt_now);
        end
    end
end

reshapedRMSE = reshape(RMSE_avg, [5, 2]);
T1 = array2table(reshapedRMSE);
reshapedRMSEpcnt = reshape(RMSEpcnt_avg, [5, 2]);
T2 = array2table(reshapedRMSEpcnt);
reshapedR2 = reshape(R2_avg, [5, 2]);
T3 = array2table(reshapedR2);
reshapedCorr = reshape(corrP_avg, [5, 2]);
T4 = array2table(reshapedCorr);

%% Plots
if strcmp(nnFile, 'Pu1_level_all_network') || strcmp(nnFile, 'Pu1_all_all_network')
    i1 = [Or3_indLG Or3_indDF Pu1_indLG Ye3_indDF];
    k1 = [3 3 5 4];    
elseif strcmp(nnFile, 'Or3_level_all_network') || strcmp(nnFile, 'Or3_all_all_network') || strcmp(nnFile, 'BL4_level_all_network') || strcmp(nnFile, 'BL4_all_all_network')
    i1 = [Or3_indLG Or3_indDF BL4_indLG BL4_indDF];
    k1 = [3 3 2 2];
elseif strcmp(nnFile, 'Bl3_level_all_network') || strcmp(nnFile, 'Bl3_all_all_network')
    i1 = [Or3_indLG Or3_indDF Bl3_indLG Bl3_indDF];
    k1 = [3 3 1 1];
end

l1 = [1 2 1 2]; 

figure
for i = 1: length(i1)
    subplot(length(i1),1,i)
    plot(time(i1(i),k1(i),l1(i)).t, force(i1(i),1,k1(i),l1(i)).meas, 'k', 'LineWidth', 2)
    hold on
    for j = 1:length(model_names)
        plot(time(i1(i),k1(i),l1(i)).t, force(i1(i),j,k1(i),l1(i)).est*f_max(i1(i),j,k1(i),l1(i)), "LineWidth", 2)
    end
    xlim([10 12])
end
legend({'Measurements', model_names{:}})

%% Check if force-length and force-velocity relationships are followed
model_params = {nnFile,'normal'};
model_names = {'NN', 'Hill'};
act = 0.2:0.2:1;
for i = 1:length(act)
    l_ce1 = 0.6:0.01:1.4;
    l_ce1 = l_ce1';
    for j = 1:length(model_names)
            flce(:,j,i) = getForce(l_ce1, zeros(size(l_ce1)), act(i)+zeros(size(l_ce1)), model_names{j}, model_params{j}, muscle_names{1});
    end
    
    v_ce1 = -10:0.1:10;
    v_ce1 = v_ce1';
    for j = 1:length(model_names)
            gvce(:,j,i) = getForce(ones(size(v_ce1)), v_ce1, act(i)+zeros(size(v_ce1)), model_names{j}, model_params{j}, muscle_names{1});
    end
end

colors = {'#f3e79b', '#fac484', '#eb7f86', '#ce6693', '#5c53a5'};
styles = {'-', '--'};
width = 2;

figure
subplot(1,2,1)
hold on
for i = 1:length(act)
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
for i = 1:length(act)
    for j = 1:length(model_names)
        plot(v_ce1, gvce(:,j,i), "Color", colors{i}, "LineStyle", styles{j}, "LineWidth", width)
        txt=[mdls{j}, ', a=' num2str(act(i))];
        legendEntries{plotCounter} = txt;
        plotCounter = plotCounter + 1;
    end
end
legend(legendEntries)
xlabel('Normalized fibre velocity')
title('Force-Velocity Relationship')


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