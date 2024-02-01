%The following code is used to train a NN model on the LG data of the Bl3
%bird. The trials that were used for training were trial "r01" and "r12".
%using the data from these  trials the training (the following code) was
%ran 10 times to get the best NN model with the lowest error, out of these
%10 iterations.  

clear all
close all
clc

%% Check path
addpath(genpath('C:\Users\annek\Documents\MATLAB\MuscleNN\musclenn\PaperCode')) %Change this to the folder where you downloaded this code

%% Settings
bird_name = 'Bl3'; 
muscle_name = 'LG'; %
trial_name = 'Bl3d2_r01_1p8_Lev_Cal'; %'Bl3d2_r12_4p5_7cm_Cal';%

%Write here the location of the data
folder_base = 'C:\Users\annek\Documents\MATLAB\MuscleModel\MuscleData\Guinnea Fowls\';

%% Load data
bird_data = readmatrix([folder_base 'MuscleMorphologyData']);

musvar = getMuscleVariables(bird_data, bird_name, muscle_name);
[time, l_ce, v_ce, EMG, Force, h] = loadDataFile(bird_name, muscle_name, trial_name,folder_base); 

%Normalize force to maximum isometric force
tendonForce =  Force/musvar.f_max; 

 % Use 80% for training
ind_use = round(length(l_ce)*0.8);

%% Network normalization to have 0 mean and std of 1 using only training data
%find the mean and stDev of the variables of the training and validation data (80-20) partitions
data_mean.EMG = mean(EMG(1:ind_use));
data_std.EMG  = std(EMG(1:ind_use));

data_mean.force = mean(tendonForce(1:ind_use));
data_std.force  = std(tendonForce(1:ind_use));

EMG_norm = doDataNormalization(EMG, data_mean.EMG,data_std.EMG);
force_norm = doDataNormalization(tendonForce, data_mean.force,data_std.force);

%% Define training and testing dataset

%Inputs
Xtrain = [EMG_norm(1:ind_use)];
%Outputs
Ytrain = force_norm(1:ind_use);

%Inputs
Xtest = [EMG_norm(ind_use+1:end)];
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

save('network_r01_emg.mat', 'NN')
