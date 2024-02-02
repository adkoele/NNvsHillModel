clear all
close all
clc

%% Add to path
%Please make sure that you use the correct \ or / for your operating system
addpath(genpath('C:\Users\annek\Documents\GitHub\NNvsHillModel\')) %Change this to the folder where you downloaded this code

%% Settings

bird_names = {'Bl3', 'BL4', 'Or3', 'Ye3', 'Pu1'}; %Birds in order of 1 to 5
muscle_names = {'LG', 'DF'};

%Write the location of the data and networks/model parameters here. Please
%make sure that you use the correct \ or / for your operating system
folder_data = 'C:\Users\annek\Desktop\GuineaFowlTest\GuineaFowlDataNNModelTraining\Guinnea Fowls Share\'; %Change this to the folder where the guinnea fowl data is available
folder_networsparams = 'C:\Users\annek\Desktop\GuineaFowlTest\NetworksAndParameterFiles\NetworksAndParameterFiles\'; %Change this to the folder where the networks and model parameter sets are stored

bird_data = readmatrix([folder_data 'MuscleMorphologyData']);
warning('The columns of the xlsx file are hard coded, please make sure that they match your version of MuscleMorphologyData.xlsx in getMuscleParameters.m')

model_params = {'network_large'};% %use 'normal' to compare to original Hill parameters and 'optresult_150523_r1.mat', 'optresult_150523_r12.mat' to compare to optimized Hill models
model_names = {'NN'};

for k = 1:length(bird_names)
    trial_names = dir([folder_data filesep bird_names{k}]);
    % Find some indices for plotting
    if strcmpi(bird_names{k}, 'or3')
        clear trial_name
        for iT = 1:length(trial_names)
            trial_name{iT} = trial_names(iT).name;
        end
        Or3_ind = find(~contains(trial_name, 'r08') == 0);
    elseif strcmpi(bird_names{k}, 'ye3')
        clear trial_name
        trial_names = dir([folder_data filesep bird_names{k}]); 
        for iT = 1:length(trial_names)
            trial_name{iT} = trial_names(iT).name;
        end
        Ye3_ind = find(~contains(trial_name, 'd2_t11') == 0);
    elseif strcmpi(bird_names{k}, 'pu1')
        clear trial_name
        trial_names = dir([folder_data filesep bird_names{k}]); 
        for iT = 1:length(trial_names)
            trial_name{iT} = trial_names(iT).name;
        end
        Pu1_ind = find(~contains(trial_name, 'r08_3p8') == 0);
    end

    for l = 1:length(muscle_names)
        for i = 3:length(trial_names)
            if ~strcmp(trial_names(i).name(end-2:end), 'mat')
                continue
            end 
            if strcmpi(bird_names{k}, 'ye3') && strcmpi(muscle_names{l}, 'df') %do not use day 1 for DF of ye3
                if strcmpi(trial_names(i).name(1:5), 'Ye3d1')
                    continue
                end
            end
            if strcmpi(bird_names{k}, 'ye3') && strcmpi(muscle_names{l}, 'lg') %do not use LG of ye3
                continue
            end
            if strcmpi(bird_names{k}, 'pu1') && strcmpi(muscle_names{l}, 'df') %do not use DF of pu1
                continue
            end
            %% Get forces, calculate RMSEs and correlations
            [time(i,k,l).t, data(i,k,l).l_ce, data(i,k,l).v_ce, data(i,k,l).EMG, force(i,1,k,l).meas] = loadDataFile(bird_names{k}, muscle_names{l}, trial_names(i).name,folder_data); 
            for j = 1:length(model_names)
                [force(i,j,k,l).est,f_max(i,j,k,l)] = getForce(folder_networsparams, data(i,k,l).l_ce, data(i,k,l).v_ce, data(i,k,l).EMG, model_names{j}, model_params{j}, muscle_names{l}, bird_names{k}, bird_data);
        
                RMSE(i,j,k,l) = rmse(force(i,j,k,l).est(:), force(i,1,k,l).meas(:)/f_max(i,j,k,l));
                corrP(i,j,k,l) = corr(force(i,j,k,l).est(:), force(i,1,k,l).meas(:)/f_max(i,j,k,l));
                [r2(i,j,k,l), RMSE1(i,j,k,l)] = rsquare(force(i,1,k,l).meas(:)/f_max(i,j,k,l),force(i,j,k,l).est(:));
            end
        end
    end
end

%% Calculate mean RMSE and correlations, not including training data
for j = 1:length(model_names)
    for k = 1:length(bird_names)
        for l = 1:length(muscle_names)
%             ind_table = (k-1)*length(muscle_names)+l;
            RMSE_now = RMSE(:,j,k,l);
            corrP_now = corrP(:, j,k,l);
            RMSE1_now = RMSE1(:,j,k,l);
            R2_now = r2(:, j,k,l);
            

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
        end
    end
end

%% Plots
trial_names = {'Or3d1_r08_3p5_7cm_Cal', 'Or3d1_r08_3p5_7cm_Cal', 'Pu1d1_r08_3p8_7cm_Cal', 'ye3d2_t11_3p5_7cm_Cal'};

i1 = [Or3_ind Or3_ind Pu1_ind Ye3_ind];
k1 = [3 3 5 4];
l1 = [1 2 1 2]; 

figure
for i = 1: length(trial_names)
    subplot(length(trial_names),1,i)
    plot(time(i1(i),k1(i),l1(i)).t, force(i1(i),1,k1(i),l1(i)).meas, 'k')
    hold on
    for j = 1:length(model_names)
        plot(time(i1(i),k1(i),l1(i)).t, force(i1(i),j,k1(i),l1(i)).est*f_max(i1(i),j,k1(i),l1(i)))
    end
    xlim([10 12])
end
legend({'Measurements', model_names{:}})

%% Check if force-length and force-velocity relationships are followed
model_params = {'network_large', 'normal'};
model_names = {'NN', 'Hill'};
act = 0.2:0.2:1;
for i = 1:length(act)
    l_ce1 = 0.6:0.01:1.4;
    l_ce1 = l_ce1';
    for j = 1:length(model_names)
            flce(:,j,i) = getForce(folder_networsparams, l_ce1, zeros(size(l_ce1)), act(i)+zeros(size(l_ce1)), model_names{j}, model_params{j}, muscle_names{1});
    end
    
    v_ce1 = -10:0.1:10;
    v_ce1 = v_ce1';
    for j = 1:length(model_names)
            gvce(:,j,i) = getForce(folder_networsparams, ones(size(v_ce1)), v_ce1, act(i)+zeros(size(v_ce1)), model_names{j}, model_params{j}, muscle_names{1});
    end
end

figure
subplot(1,2,1)
hold on
for i = 1:length(act)
    for j = 1:length(model_names)
        plot(l_ce1, flce(:,j,i))
    end
end
xlabel('Normalized fibre length')
ylabel('Normalized muscle force')
title('Force-Length Relationship')

subplot(1,2,2)
hold on
for i = 1:length(act)
    for j = 1:length(model_names)
        plot(v_ce1, gvce(:,j,i))
    end
end
legend('Hill-type model', 'Neural Network')
xlabel('Normalized fibre velocity')
title('Force-Velocity Relationship')