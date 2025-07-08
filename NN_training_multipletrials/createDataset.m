%The following code creates a training and validation dataset based on
%user inputs which define the inclusion criteria for the trials
%(i.e. trial speed and obstacle height). The resulting dataset will be
%used for training a neural network.

clear all
close all
clc

%Note that running this file takes very long. Therefore, it is recommended
%to do it only once and use the .mat file after.

%% Check path
filenamePath=mfilename('fullpath');
filePath =[fileparts(filenamePath) filesep '..' filesep];
addpath(genpath(filePath)) 

%% Settings
bird_names = {}; %Will be set according to user input 
obstacleHeight = "";
speed = "";
muscle_names = {'LG'}; %The muscle(s) on which to train the neural network

%Write here the location of the data
folder_base = '/path/to/NN_Data/';
bird_data = readmatrix([folder_base 'MuscleMorphologyData']);

lcetrain_all = []; vcetrain_all = []; EMGtrain_all = []; Forcetrain_all = [];

lceval_all = []; vceval_all = []; EMGval_all = []; Forceval_all = [];

[lce_train, lce_val, vce_train, vce_val, EMG_train, EMG_val, Force_train, Force_val, t_train, t_val] = deal(cell(20, length(bird_names)));

tTrain_all = [];
tVal_all = [];
filenames_all = {};

%% User input requested for the generation of the training and validation dataset
% User options initialized with erroneous values
excludeBird=4;
excludeObstacleHeight=2;
speedOption=2;

while ~ismember(excludeBird, [1 2 3 5])
    excludeBird = input('Please enter: 1 to exclude Bird 1/Bl3, 2 to exclude Bird 2/Bl4, 3 to exclude Bird 3/Or3 and 5 to exclude Bird 5/Pu1. ');
    if  ~ismember(excludeBird, [1 2 3 5])
        disp('Invalid input. Please enter either 1, 2, 3 or 5.')
    end
end

while ~ismember(excludeObstacleHeight, [0 1])
    excludeObstacleHeight = input('Please enter: 0 to use trials with any and all obstacle heights or 1 to use only level trials. ');
    if  ~ismember(excludeObstacleHeight, [0 1])
        disp('Invalid input. Please enter either 0 or 1.')
    end
end

while ~ismember(speedOption, [0 1])
    speedOption = input('Please enter: 0 to use trials with any and all speeds or 1 to only use trials with the lowest speed (1.8m/s). ');
    if  ~ismember(speedOption, [0 1])
        disp('Invalid input. Please enter either 0 or 1.')
    end
end

% Setting the names of the birds used for training and of the excluded 
% bird, according to the previous user input
switch excludeBird
    case {1}
        bird_names={'BL4', 'Or3', 'Pu1'};
        excluded='Bl3';
    case {2}
        bird_names={'Bl3', 'Or3', 'Pu1'};
        excluded='BL4';
    case {3}
        bird_names={'Bl3', 'BL4', 'Pu1'};
        excluded='Or3';
    case {5}
        bird_names={'Bl3', 'BL4', 'Or3'};
        excluded='Pu1';
    otherwise
        error('Value Error: Expected value 1, 2, 3 or 5');
end

% Setting the obstacle height, according to the previous user input
switch excludeObstacleHeight
    case {0}
        obstacleHeight='all';
    case {1}
        obstacleHeight='level';
    otherwise
        error('Value Error: Expected value 0 or 1');
end

% Setting the speed, according to the previous user input. Higher speeds
% than 1.8 are not implemented
switch speedOption
    case {0}
        speed='all';
    case {1}
        speed='1p8';
    otherwise
        error('Value Error: Expected value 0 or 1');
end

% Storing the user-defined settings in the object that will also store the
% generated network. These saved settings/parameters are later on
% needed and used in the analysis code.
paramStruct.birds = bird_names;
paramStruct.excludedBird = excluded;
paramStruct.obstacleHeight = obstacleHeight;
paramStruct.speed = speed;

% Selecting 1 trial, from each of the 3 birds which will be used for training
% the NN, so that it is reserved for testing. Different trials are selected
% depending on the user-specified dataset characteristics (i.e. level
% trials or trials with obstacles, etc.)
if strcmpi(obstacleHeight,'all') && strcmpi(speed, 'all')
    substrings = {};
    for ind = 1:length(bird_names)
        if strcmpi(bird_names{ind}, 'Bl3')
            trialsForTesting{ind}='Bl3d2_r07_4p5_5cm_Cal.mat';
        elseif strcmpi(bird_names{ind}, 'BL4')
            trialsForTesting{ind}='BL4d2_r02_3p8_Lev_Cal.mat';
        elseif strcmpi(bird_names{ind}, 'Or3')
            trialsForTesting{ind}='Or3d1_r05_1p8_7cm_Cal.mat';
        elseif strcmpi(bird_names{ind}, 'Pu1')
            trialsForTesting{ind}='Pu1d1_r05_1p8_5cm_Cal.mat';
        end
    end
    paramStruct.testTrials=trialsForTesting;
elseif strcmpi(obstacleHeight,'level') && strcmpi(speed, 'all')
    substrings = {'5cm', '7cm'};
    for ind = 1:length(bird_names)
        if strcmpi(bird_names{ind}, 'Bl3')
            trialsForTesting{ind}='Bl3d2_r03_4p5_Lev_Cal.mat';
        elseif strcmpi(bird_names{ind}, 'BL4')
            trialsForTesting{ind}='BL4d2_r02_3p8_Lev_Cal.mat';
        elseif strcmpi(bird_names{ind}, 'Or3')
            trialsForTesting{ind}='Or3d1_r03_1p8_Lev_Cal.mat';
        elseif strcmpi(bird_names{ind}, 'Pu1')
            trialsForTesting{ind}='Pu1d1_r01_1p8_Lev_Cal.mat';
        end
    end
    paramStruct.testTrials=trialsForTesting;
elseif strcmpi(obstacleHeight,'all') && strcmpi(speed, '1p8')
    substrings = {'3p8', '4p5', '3p0', '3p5'};
    for ind = 1:length(bird_names)
        if strcmpi(bird_names{ind}, 'Bl3')
            trialsForTesting{ind}='Bl3d2_r01_1p8_Lev_Cal.mat';
        elseif strcmpi(bird_names{ind}, 'BL4')
            trialsForTesting{ind}='BL4d2_r10_1p8_5cm_Cal.mat';
        elseif strcmpi(bird_names{ind}, 'Or3')
            trialsForTesting{ind}='Or3d1_r05_1p8_7cm_Cal.mat';
        elseif strcmpi(bird_names{ind}, 'Pu1')
            trialsForTesting{ind}='Pu1d1_r05_1p8_5cm_Cal.mat';
        end
    end
    paramStruct.testTrials=trialsForTesting;
elseif strcmpi(obstacleHeight,'level') && strcmpi(speed, '1p8')
    substrings = {'5cm', '7cm', '3p8', '4p5', '3p0', '3p5'};
    trialsForTesting={};
    paramStruct.testTrials="No test trials reserved. Training data too few.";
end

for iFolder=1:length(bird_names)
    directoryStr= [folder_base bird_names{iFolder} filesep];
    directoryInstance=dir(directoryStr);
    filenames={directoryInstance.name};
        
         
   for iMus = 1:length(muscle_names)
        muscleParams = getMuscleParameters(bird_data, bird_names{iFolder}, muscle_names{iMus});
        f_max = muscleParams.f_max;

        
        % Setting the final list of trial names that will be used for
        % training the NN
        if isempty(substrings) && ~isempty(trialsForTesting)
            filenames = setdiff(filenames, trialsForTesting, 'stable');
        elseif ~isempty(substrings) && ~isempty(trialsForTesting)
            %remove trials if they contain steps (they are not level)
            mask = contains(filenames, substrings);
            filenames(mask) = [];

            %remove 3 level trials (to use later for testing the NN)
            filenames = setdiff(filenames, trialsForTesting, 'stable');
        elseif ~isempty(substrings) && isempty(trialsForTesting)
            %remove trials if they contain steps (they are not level) and
            %higher speeds
            mask = contains(filenames, substrings);
            filenames(mask) = [];
        else
            error('Recheck the contents of variables substrings and trialsForTesting.');
        end

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
              
saveAs = [filePath 'data' filesep excluded, '_' , obstacleHeight, '_' ,speed, '_dataTrainValidation.mat'];
save(saveAs, 'filenames_all', 'lcetrain_all', 'vcetrain_all', 'EMGtrain_all', 'Forcetrain_all', 'lceval_all', 'vceval_all', 'EMGval_all', 'Forceval_all', 'tTrain_all', 'tVal_all', "paramStruct")

%Function to split training and validation data
function [out_train, out_val] = separateTestValData(in, sep)
    if nargin == 1
        sep = 0.8; %use 80% of data for training
    end
    ind_use = round(length(in)*sep);
    out_train = in(1:ind_use,:);
    out_val  = in(ind_use+1:end,:);
end