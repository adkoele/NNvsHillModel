clear all
close all
clc

%% Add to path
%Please make sure that you use the correct \ or / for your operating system
addpath(genpath('/path/to/code/')) %Change this to the folder where you downloaded this code

filenamePath=mfilename('fullpath');
filePath =[fileparts(filenamePath) filesep '..' filesep];
addpath(genpath(filePath)) 
    


%% Settings
bird_name = 'Bl3'; %
muscle_name = 'LG'; %'DF'; %
trial_name = 'Bl3d2_r12_4p5_7cm_Cal'; %'Bl3d2_r01_1p8_Lev_Cal' 'Bl3d2_r12_4p5_7cm_Cal';

%Write the location of the data and networks/model parameters here. Please
%make sure that you use the correct \ or / for your operating system
folder_data = '/path/to/NN_Data/'; %Change this to the folder where the guinnea fowl data is available

bird_data = readmatrix([folder_data 'MuscleMorphologyData']);
warning('The columns of the xlsx file are hard coded, please make sure that they match your version of MuscleMorphologyData.xlsx in getMuscleParameters.m')

%% Load data
musvar = getMuscleParameters(bird_data, bird_name, muscle_name);
[time, l_ce, v_ce, EMG, Force, h] = loadDataFile(bird_name, muscle_name, trial_name, folder_data); 

%% Generate standard muscle model (parameters to be optimized later)
modelvar.v_max = 10;
modelvar.PEEslack = 1.2; 
modelvar.gmax = 1.5; 
modelvar.kPEE = 1/musvar.l_opt^2;
modelvar.Arel = 0.25; 
modelvar.W = 0.4; %

%dependent parameters         
modelvar.c3 = modelvar.v_max*modelvar.Arel*(modelvar.gmax - 1.)/(modelvar.Arel + 1);

%% Perform optimization
%First 80% for optimization
ind_use = round(length(l_ce)*0.8);

optvar = cmaes(modelvar, musvar, l_ce(1:ind_use), v_ce(1:ind_use), EMG(1:ind_use), Force(1:ind_use));

dateString = datestr(now, 'ddmmyyyy_HHMM');
    
saveAs =  [filePath 'Hill_optimization/' dateString '_' trial_name '_baselineHillModel.mat']; %naming convension: nameOfExcludedBird_obstacleHeight_speed_...
save(saveAs, 'optvar');

disp('Saved baseline Hill model');