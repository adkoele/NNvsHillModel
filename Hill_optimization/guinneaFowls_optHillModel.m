clear all
close all
clc

%% Check path
addpath(genpath('C:\Users\annek\Documents\MATLAB\MuscleNN\musclenn\PaperCode')) %Change this to the folder where you downloaded this code

%% Settings
bird_name = 'Bl3'; %
muscle_name = 'LG'; %'DF'; %
trial_name = 'Bl3d2_r01_1p8_Lev_Cal';% 'Bl3d2_r12_4p5_7cm_Cal';

%Write here the location of the data
folder_base = 'C:\Users\annek\Documents\MATLAB\MuscleModel\MuscleData\Guinnea Fowls\'; %Change this folder to where the data is located

%% Load data
musvar = getMuscleVariables(folder_base, bird_name, muscle_name);
[time, l_ce, v_ce, EMG, Force, h] = loadDataFile(bird_name, muscle_name, trial_name, folder_base); 

%% Generate standard muscle model (parameters to be optimized later)
modelvar.v_max = 10;
modelvar.PEEslack = 1.2;
modelvar.gmax = 1.5;
modelvar.kPEE = 1/musvar.l_opt^2;
modelvar.Arel = 0.25;
modelvar.W = 0.4;

%dependent parameters         
modelvar.c3 = modelvar.v_max*modelvar.Arel*(modelvar.gmax - 1.)/(modelvar.Arel + 1);

%% Perform optimization
%First 80% for optimization
ind_use = round(length(l_ce)*0.8);

optvar = cmaes(modelvar, musvar, l_ce(1:ind_use), v_ce(1:ind_use), EMG(1:ind_use), Force(1:ind_use));