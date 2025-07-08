% A high-pass Butterworth filter (30 Hz) was applied to remove low frequency noise, the signal was
% rectified, and a low-pass Butterworth filter (6 Hz) was applied to generate a smoothed EMG signal
% using a moving average zero-phase digital filter with a window of 120 ms.
function filteredEMG = filtEMG(dta, muscleIndx, frq)

    EMG2 = dta.mV_EMG(:,muscleIndx);
    [B,A] = butter(3,30/(frq/2),'high');
    EMG3 = filtfilt(B,A,EMG2); %high-pass filter
    EMG4 = abs(EMG3); %rectification
    [B,A] = butter(3,6/(frq/2),'low');
    filteredEMG = max(filtfilt(B,A,EMG4),0); %low-pass filter
