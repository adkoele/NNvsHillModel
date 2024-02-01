clear all
close all
clc

%Note that running this file takes very long. Therefore, it is recommended
%to do it only once and use the .mat file after.

%% Check path
addpath(genpath('C:\Users\annek\Documents\MATLAB\MuscleNN\musclenn\PaperCode')) %Change this to the folder where you downloaded this code

%% Settings
bird_names = {'Bl3', 'BL4', 'Or3'}; 
muscle_names = {'LG'};%, 'DF'}; %
trial_names = 'Bl3d2_r01_1p8_Lev_Cal'; %'Bl3d2_r12_4p5_7cm_Cal';%

%Write here the location of the data
folder_base = 'C:\Users\annek\Documents\MATLAB\MuscleModel\MuscleData\Guinnea Fowls\';
bird_data = readmatrix([folder_base 'MuscleMorphologyData']);

lcetrain_all = []; vcetrain_all = []; EMGtrain_all = []; Forcetrain_all = [];

lceval_all = []; vceval_all = []; EMGval_all = []; Forceval_all = [];

[lce_train, lce_val, vce_train, vce_val, EMG_train, EMG_val, Force_train, Force_val, t_train, t_val] = deal(cell(20, length(bird_names)));

tTrain_all = [];
tVal_all = [];
filenames_all = {};

for iFolder=1:length(bird_names)
    directoryStr= [folder_base bird_names{iFolder} filesep];
    directoryInstance=dir(directoryStr);
    filenames={directoryInstance.name};
        
   maximumStress=0.36;
     
   for iMus = 1:length(muscle_names)
        ind_row = findBird(bird_names{iFolder});
        %birdData
        if strcmpi(muscle_names{iMus}, 'lg')
            PCSA = bird_data(ind_row,8); %converted to m
        elseif strcmpi(muscle_names{iMus}, 'df')
            PCSA = bird_data(ind_row,14);
        else
            error('Incorrect muscle name')
        end
        
        f_max=PCSA*maximumStress; %multiply the PCSA of tendon G with maximum stress
        
        iFile = 1;
        for i = 3:length(filenames)
            usedFile=filenames{i};
            dataFile=directoryStr+""+usedFile;
        
            if strcmpi(usedFile(end-2:end), 'mat')
                load(dataFile);
               
                %check filenames
                [time, l_ce, v_ce, EMG, Forcetrain_all, h] = loadDataFile(bird_names{iFolder}, muscle_names{iMus}, filenames{i},folder_base); 
                
                %Normalize force to maximum isometric force
                Force_norm =  Forcetrain_all/f_max; 
                %------------------------------------------------
                % DATA PARTITIONING    
                [lce_train{iFile, iFolder}, lce_val{iFile, iFolder}] = separateTestValData(l_ce);
                [vce_train{iFile, iFolder}, vce_val{iFile, iFolder}] = separateTestValData(v_ce);
                [EMG_train{iFile, iFolder}, EMG_val{iFile, iFolder}] = separateTestValData(EMG);
                [Force_train{iFile, iFolder}, Force_val{iFile, iFolder}] = separateTestValData(Force_norm);
                [t_train{iFile, iFolder}, t_val{iFile, iFolder}] = separateTestValData(time);
                
                %-----------------------------------------------

                filenames_all{iFolder, iFile} = usedFile;

                iFile = iFile+1;
            end
        end
   end
end 

lcetrain_all = vertcat(lce_train{:});
vcetrain_all = vertcat(vce_train{:});
EMGtrain_all = vertcat(EMG_train{:});
Forcetrain_all = vertcat(Force_train{:});
tTrain_all = vertcat(t_train{:});

lceval_all = vertcat(lce_val{:});
vceval_all = vertcat(vce_val{:});
EMGval_all = vertcat(EMG_val{:});
Forceval_all = vertcat(Force_val{:});
tVal_all = vertcat(t_val{:});
              
save('data_train_validation.mat', 'filenames_all', 'lcetrain_all', 'vcetrain_all', 'EMGtrain_all', 'Forcetrain_all', 'lceval_all', 'vceval_all', 'EMGval_all', 'Forceval_all', 'tTrain_all', 'tVal_all')

function [out_train, out_val] = separateTestValData(in, sep)
    if nargin == 1
        sep = 0.8; %use 80% of data for training
    end
    ind_use = round(length(in)*sep);
    out_train = in(1:ind_use,:);
    out_val  = in(ind_use+1:end,:);
end