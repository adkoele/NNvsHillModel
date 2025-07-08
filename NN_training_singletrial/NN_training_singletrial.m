%The following code is used to train a NN model on the LG data of the Bl3
%bird. The trials that were used for training were trial "r01" and "r12".
%using the data from these  trials the training (the following code) was
%ran 10 times to get the best NN model with the lowest error, out of these
%10 iterations.  

clear all
close all
clc

%% Check path
addpath(genpath('/path/to/code/')) %Change this to the folder where you downloaded this code

%% Settings
bird_name = 'Bl3'; 
muscle_name = 'LG'; 
trial_name = {'Bl3d2_r12_4p5_7cm_Cal'}; % Options: 'Bl3d2_r01_1p8_Lev_Cal'; 'Bl3d2_r12_4p5_7cm_Cal';%

dateString = datestr(date, 'dd_mm_yyyy');

%Write here the location of the data
folder_base = '/path/to/NN_Data/';

%% Load data
bird_data = readmatrix([folder_base 'MuscleMorphologyData']);

for trialI=1:length(trial_name)

    musvar = getMuscleParameters(bird_data, bird_name, muscle_name);
    [time, l_ce, v_ce, EMG, Force, h] = loadDataFile(bird_name, muscle_name, trial_name{trialI}, folder_base); 
    
    %Normalize force to maximum isometric force
    tendonForce =  Force/musvar.f_max; 
    
     % Use 80% for training
    ind_use = round(length(l_ce)*0.8);
    
    %% Network normalization to have 0 mean and std of 1 using only training data
    %find the mean and stDev of the variables of the training and validation data (80-20) partitions
    data_mean.lce = mean(l_ce(1:ind_use));
    data_std.lce  = std(l_ce(1:ind_use));
    
    data_mean.vce = mean(v_ce(1:ind_use));
    data_std.vce  = std(v_ce(1:ind_use));
    
    data_mean.EMG = mean(EMG(1:ind_use));
    data_std.EMG  = std(EMG(1:ind_use));
    
    data_mean.force = mean(tendonForce(1:ind_use));
    data_std.force  = std(tendonForce(1:ind_use));
    
    
    lce_norm = doDataNormalization(l_ce, data_mean.lce,data_std.lce);
    vce_norm = doDataNormalization(v_ce, data_mean.vce,data_std.vce);
    EMG_norm = EMG; 
    force_norm = doDataNormalization(tendonForce, data_mean.force,data_std.force);
    
    %% Define training and testing dataset
    
    %Inputs
    Xtrain = [lce_norm(1:ind_use), vce_norm(1:ind_use), EMG_norm(1:ind_use)];
    %Outputs
    Ytrain = force_norm(1:ind_use);
    
    %Inputs
    Xtest = [lce_norm(ind_use+1:end), vce_norm(ind_use+1:end), EMG_norm(ind_use+1:end)];
    %Outputs
    Ytest = force_norm(ind_use+1:end);
    
    %% Training with hyperparameter optimization
    
    rng("default") % For reproducibility
    
    for i = 1:10 %Repeat ten times
        res(i).Mdl = fitrnet(Xtrain, Ytrain,'OptimizeHyperparameters', 'auto', "HyperparameterOptimizationOptions", struct("MaxObjectiveEvaluations", 60),"ValidationData", {Xtest, Ytest});
        loss_val(i) = res(i).Mdl.HyperparameterOptimizationResults.MinObjective;
    end
    [~,iMin] = min(loss_val); %Find lowest objective
    Mdl = res(iMin).Mdl;
    
    testMSE = loss(Mdl,Xtest,Ytest);
    trainMSE = loss(Mdl,Xtrain,Ytrain);
    
    %% Plot training data
    
    iteration = Mdl.TrainingHistory.Iteration;
    trainLosses = Mdl.TrainingHistory.TrainingLoss;
    valLosses = Mdl.TrainingHistory.ValidationLoss;
    plot(iteration,trainLosses,iteration,valLosses)
    legend(["Training","Validation"])
    xlabel("Iteration")
    ylabel("Mean Squared Error")
    %------------------------------------
    
    %% Save model and normalization info
    NN.Mdl = Mdl;
    NN.data_mean = data_mean;
    NN.data_std = data_std;
    
    if contains(trial_name{trialI}, 'r01')
        netwName=['network_r01_' dateString];
    elseif contains(trial_name{trialI}, 'r12')
        netwName=['network_r12_' dateString];
    end
    
    save([netwName '.mat'], 'NN', 'res', 'loss_val')

end
