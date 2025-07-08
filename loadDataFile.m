function [time, l_ce, v_ce, EMG_out, Force, h] = loadDataFile(bird_name, muscle_name, trialname, folder_base)

if strcmpi(bird_name,'pu1')
    folder = [folder_base 'Pu1' filesep];
elseif strcmpi(bird_name,'ye3')
    folder = [folder_base 'Ye3' filesep];
elseif strcmpi(bird_name,'or3')
    folder = [folder_base 'Or3' filesep];
elseif strcmpi(bird_name,'bl3')
    folder = [folder_base 'Bl3' filesep];
elseif strcmpi(bird_name,'bl4')
    folder = [folder_base 'BL4' filesep];
end

load([folder trialname]);
time = data.time - data.time(1,1);
h = mean(diff(time)); %timestep of data
freq = 1/h;
if strcmpi(muscle_name, 'lg')
    ind_mus = 1;
elseif strcmpi(muscle_name, 'df')
    ind_mus = 2;
else
    error('Incorrect muscle name')
end
l_ce = data.muscleLength_FL(:,ind_mus);
v_ce_raw = data.velocity_FL(:,ind_mus);

% filter v_ce
[B,A] = butter(3,6/(freq/2));
v_ce = filtfilt(B,A,v_ce_raw);

Force = data.tendonForceN(:,ind_mus);

% A high-pass Butterworth filter (30 Hz) was applied to remove low frequency noise, the signal was
% rectified, and a low-pass Butterworth filter (6 Hz) was applied to generate a smoothed EMG signal
% using a moving average zero-phase digital filter with a window of 120 ms.
EMG = filtEMG(data, ind_mus, freq);

% Check if average maximum Activation file already exists, if not generate
% it
filePath = fileparts(mfilename('fullpath'));
avgMaxActFilePath=[filePath filesep 'avgMaxActivation.mat'];

if exist(avgMaxActFilePath, 'file')==2
    disp('File with average maximum activation, from max speed trials, for each bird and muscle exists');  
else
    disp('File with average maximum activation, from max speed trials, for each bird and muscle does not exist. It will be generated.');
    generateAvgMaxEMGForEachBird();   
end

% Load file with average max activation and choose correct value based on
% bird and muscle names
avgMaxActFile=load([filePath filesep 'avgMaxActivation.mat']);

if strcmpi(bird_name, 'Bl3') && ind_mus==1
    avgMaxActIndx=1;
elseif strcmpi(bird_name, 'Bl3') && ind_mus==2
    avgMaxActIndx=5;
elseif strcmpi(bird_name, 'BL4') && ind_mus==1
    avgMaxActIndx=2;
elseif strcmpi(bird_name, 'BL4') && ind_mus==2
    avgMaxActIndx=6;
elseif strcmpi(bird_name, 'Or3') && ind_mus==1
    avgMaxActIndx=3;
elseif strcmpi(bird_name, 'Or3') && ind_mus==2
    avgMaxActIndx=7;
elseif strcmpi(bird_name, 'Ye3') && ind_mus==2
    avgMaxActIndx=8;
elseif strcmpi(bird_name, 'Pu1') && ind_mus==1
    avgMaxActIndx=4;
end

avgMaxEMG=avgMaxActFile.avgMaxActPerBirdArray{avgMaxActIndx,3};



%=============================================================
% Average max normalization, that uses the average of all the max
% EMG values across highest-speed trials of
% that bird's muscle
%=============================================================
EMG1 = (EMG/avgMaxEMG); 
%=============================================================


%include time delay in EMG
td = 23.6/1000;
ind_delay = fix(td/h);

EMG_out = zeros(size(EMG1));
EMG_out(ind_delay+1:end) = EMG1(1:end-ind_delay);
