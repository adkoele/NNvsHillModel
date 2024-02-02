clear all
close all
clc

%% Add to path
%Please make sure that you use the correct \ or / for your operating system
addpath(genpath('C:\Users\annek\Documents\GitHub\NNvsHillModel\')) %Change this to the folder where you downloaded this code

%% Settings

bird_names = {'Bl3', 'BL4'};
muscle_names = {'LG', 'DF'};

%Write the location of the data and networks/model parameters here. Please
%make sure that you use the correct \ or / for your operating system
folder_data = 'C:\Users\annek\Desktop\GuineaFowlTest\GuineaFowlDataNNModelTraining\Guinnea Fowls Share\'; %Change this to the folder where the guinnea fowl data is available
folder_networsparams = 'C:\Users\annek\Desktop\GuineaFowlTest\NetworksAndParameterFiles\NetworksAndParameterFiles\'; %Change this to the folder where the networks and model parameter sets are stored

bird_data = readmatrix([folder_data 'MuscleMorphologyData']);
warning('The columns of the xlsx file are hard coded, please make sure that they match your version of MuscleMorphologyData.xlsx in getMuscleParameters.m')

model_params = {'optresult_150523_r1.mat', 'optresult_150523_r12.mat', 'network_r01' 'network_r12'};
model_names = {'Hill', 'Hill', 'NN', 'NN'};

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

    for l = 1:length(muscle_names)
        for i = 3:length(trial_names)
            if ~strcmp(trial_names(i).name(end-2:end), 'mat')
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
            ind_table = (k-1)*length(muscle_names)+l;
            RMSE_now = RMSE(:,j,k,l);
            corrP_now = corrP(:, j,k,l);
            RMSE1_now = RMSE1(:,j,k,l);
            R2_now = r2(:, j,k,l);

            % identify and remove zeros
            ind_rem = find(abs(RMSE_now)<1e-6);
            %remove training data
            if strcmpi(bird_names{k}, 'bl3')
                cur_modname = model_params{j};
                r1_loc = strfind(cur_modname, 'r1');
                if strcmpi(cur_modname(r1_loc+2), '2')
                    ind_rem = [ind_rem; r12_ind];
                else
                    ind_rem = [ind_rem; r01_ind];
                end
            end

            RMSE_now(ind_rem) = [];
            RMSE_avg(ind_table,j) = mean(RMSE_now);

            corrP_now(ind_rem) = [];
            corrP_avg(ind_table,j) = mean(corrP_now);

            RMSE1_now(ind_rem) = [];
            RMSE1_avg(ind_table,j) = mean(RMSE1_now);

            R2_now(ind_rem) = [];
            R2_avg(ind_table,j) = mean(R2_now);
        end
    end
end


%% Plots
trial_names = {'Bl3d2_r10_3p8_7cm_Cal', 'Bl3d2_r10_3p8_7cm_Cal', 'BL4d2_r09_3p8_7cm_Cal', 'BL4d2_r09_3p8_7cm_Cal'};
i1 = [r10_ind r10_ind r09_ind r09_ind];
k1 = [1 1 2 2];
l1 = [1 2 1 2];


figure
for i = 1: length(trial_names)
    subplot(length(trial_names),1,i)
    plot(time(i1(i),k1(i),l1(i)).t, force(i1(i),1,k1(i),l1(i)).meas, 'k')
    hold on
    for j = 1:length(model_names)
        plot(time(i1(i),k1(i),l1(i)).t, force(i1(i),j,k1(i),l1(i)).est*f_max(i1(i),j,k1(i), l1(i)))
    end
    xlim([7.4 8.9])
end