%The following code is used to train a NN model on the LG data of three birds
% (Bl3, BL4, Or3). All trials were used for training. The training was
%repeated 10 times to get the best NN model with the lowest error, out of these
%10 iterations.  

clear all
close all
clc

%% Check path
addpath(genpath('C:\Users\annek\Documents\MATLAB\MuscleNN\musclenn\PaperCode')) %Change this to the folder where you downloaded this code

%% Settings
bird_names = {'Bl3', 'BL4', 'Or3'}; 
muscle_name = 'LG'; %

%Write here the location of the data
folder_base = 'C:\Users\annek\Documents\MATLAB\MuscleModel\MuscleData\Guinnea Fowls\';

%% Load data
load('data_train_validation.mat'); %If this file does not yet exist, please run createDataset.m once first

%% Network normalization to have 0 mean and std of 1 using only training data
%find the mean and stDev of the variables of the training and validation data (80-20) partitions
data_mean.lce = mean(lcetrain_all);
data_std.lce  = std(lcetrain_all);

data_mean.vce = mean(vcetrain_all);
data_std.vce  = std(vcetrain_all);

data_mean.EMG = mean(EMGtrain_all);
data_std.EMG  = std(EMGtrain_all);

data_mean.force = mean(Forcetrain_all);
data_std.force  = std(Forcetrain_all);

lcetrain_norm = doDataNormalization(lcetrain_all, data_mean.lce,data_std.lce);
vcetrain_norm = doDataNormalization(vcetrain_all, data_mean.vce,data_std.vce);
EMGtrain_norm = doDataNormalization(EMGtrain_all, data_mean.EMG,data_std.EMG);
Forcetrain_norm = doDataNormalization(Forcetrain_all, data_mean.force,data_std.force);

lceval_norm = doDataNormalization(lceval_all, data_mean.lce,data_std.lce);
vceval_norm = doDataNormalization(vceval_all, data_mean.vce,data_std.vce);
EMGval_norm = doDataNormalization(EMGval_all, data_mean.EMG,data_std.EMG);
Forceval_norm = doDataNormalization(Forceval_all, data_mean.force,data_std.force);

%% Define training and testing dataset

%Inputs
Xtrain = [lcetrain_norm vcetrain_norm EMGtrain_norm];
%Outputs
Ytrain = Forcetrain_norm;

%Inputs
Xval = [lceval_norm vceval_norm EMGval_norm Forceval_norm];
%Outputs
Yval = Forceval_norm;

%% Training with hyperparameter optimization

rng("default") % For reproducibility

for i = 1:10 %Repeat ten times
    res(i).Mdl = fitrnet(Xtrain, Ytrain,'OptimizeHyperparameters', 'auto', "HyperparameterOptimizationOptions", struct("MaxObjectiveEvaluations", 60),"ValidationData", {Xval, Yval});
    loss(i) = res(i).Mdl.HyperparameterOptimizationResults.MinObjective;
end
[~,iMin] = min(loss); %Find lowest objective
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

save('network_large.mat', 'NN')
