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
EMG2 = data.mV_EMG(:,ind_mus);
[B,A] = butter(3,30/(freq/2),'high');
EMG3 = filtfilt(B,A,EMG2); %high-pass filter
EMG4 = abs(EMG3); %rectification
[B,A] = butter(3,6/(freq/2),'low');
EMG = max(filtfilt(B,A,EMG4),0); %low-pass filter

% Load all files to find max recorded value:
list_files = dir(folder);
EMG_all = [];
for i = 3:length(list_files)
    if strcmp(list_files(i).name(end-2:end), 'mat')
        load([folder list_files(i).name]);
    else
        continue
    end
    EMG2 = data.mV_EMG(:,ind_mus);
    [B,A] = butter(3,30/(freq/2),'high');
    EMG3 = filtfilt(B,A,EMG2); %high-pass filter
    EMG4 = abs(EMG3); %rectification
    [B,A] = butter(3,6/(freq/2),'low');
    EMG_temp = max(filtfilt(B,A,EMG4),0);%filtfilt(B,A,EMG4); %low-pass filter
    EMG_all = [EMG_all;EMG_temp];
end
EMG_max = max(EMG_all);
EMG1 = EMG/EMG_max;

%include time delay in EMG
td = 23.6/1000;
ind_delay = fix(td/h);

EMG_out = zeros(size(EMG1));
EMG_out(ind_delay+1:end) = EMG1(1:end-ind_delay);
