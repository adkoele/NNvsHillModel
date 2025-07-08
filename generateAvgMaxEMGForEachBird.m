function generateAvgMaxEMGForEachBird()

%% Check path
addpath(genpath('/path/to/code/')) %Change this to the folder where you downloaded this code

%% Settings
bird_namesLG = {'Bl3', 'BL4', 'Or3', 'Pu1'};  
bird_namesDF = {'Bl3', 'BL4', 'Or3', 'Ye3'}; 
elevation = "";
speed = "";
muscle_names = {'LG', 'DF'}; 

%Write here the location of the data
folder_base = '/path/to/NN_Data/';
bird_data = xlsread([folder_base 'MuscleMorphologyData']);

filenamePath=mfilename('fullpath');
filePath=fileparts(filenamePath);


countIter=0;

for iMus = 1:length(muscle_names)
    bird_names=[];
    if strcmp(muscle_names{iMus},'LG')
        bird_names=bird_namesLG;
    elseif strcmp(muscle_names{iMus},'DF')
        bird_names=bird_namesDF;
    end
     
   for iFolder=1:length(bird_names)
        directoryStr= [folder_base bird_names{iFolder} filesep];
        directoryInstance=dir(directoryStr);
        filenames={directoryInstance.name};
        countIter=countIter+1;
             
        
        count=0;
        avgMaxActPerTrial=[];
        iFile = 1;
        for i = 3:length(filenames)
            usedFile=filenames{i};
            dataFile=directoryStr+""+usedFile;
        
            if strcmpi(usedFile(end-2:end), 'mat')
                dataStruct = load(dataFile);
                data=dataStruct.data;
                
                time = data.time - data.time(1,1);
                h = mean(diff(time)); %timestep of data
                freq = 1/h;

                
                if strcmp(bird_names{iFolder}, 'Bl3') && contains(usedFile, '4p5')
                    count=count+1;
                    EMGfiltered = filtEMG(data, iMus, freq);
                    [actAtPeak, timeAtPeak] = findpeaks(EMGfiltered, 'MinPeakHeight', 0.1);
                    avgMaxActPerTrial(count)=mean(actAtPeak);
                    
                    
                elseif strcmp(bird_names{iFolder}, 'BL4') && contains(usedFile, '4p5')
                    count=count+1;
                    EMGfiltered = filtEMG(data, iMus, freq);
                    [actAtPeak, timeAtPeak] = findpeaks(EMGfiltered, 'MinPeakHeight', 0.1);
                    avgMaxActPerTrial(count)=mean(actAtPeak);

                elseif strcmp(bird_names{iFolder}, 'Or3') && contains(usedFile, '3p5')
                    
                    count=count+1;
                    EMGfiltered = filtEMG(data, iMus, freq);
                    if iMus == 2
                        [actAtPeak, timeAtPeak] = findpeaks(EMGfiltered, 'MinPeakHeight', 0.035);
                    elseif iMus == 1
                        [actAtPeak, timeAtPeak] = findpeaks(EMGfiltered, 'MinPeakHeight', 0.1);
                    end
                    
                    avgMaxActPerTrial(count)=mean(actAtPeak);

                elseif strcmp(bird_names{iFolder}, 'Ye3') && contains(usedFile, '3p5') && contains(usedFile, 'd2') 
                    
                    count=count+1;
                    EMGfiltered = filtEMG(data, iMus, freq);
                    [actAtPeak, timeAtPeak] = findpeaks(EMGfiltered, 'MinPeakHeight', 0.1);
                    avgMaxActPerTrial(count)=mean(actAtPeak);

                elseif strcmp(bird_names{iFolder}, 'Pu1') && contains(usedFile, '3p8')
                    
                    count=count+1;
                    EMGfiltered = filtEMG(data, iMus, freq);
                    [actAtPeak, timeAtPeak] = findpeaks(EMGfiltered, 'MinPeakHeight', 0.1);
                    avgMaxActPerTrial(count)=mean(actAtPeak);

                end
               
            end            
        end
        
        avgMaxActPerBirdArray{countIter, 1} = bird_names{iFolder};
        avgMaxActPerBirdArray{countIter, 2} = muscle_names{iMus};
        avgMaxActPerBirdArray{countIter, 3} = mean(avgMaxActPerTrial);
        avgMaxActPerBird(iFolder, iMus) = mean(avgMaxActPerTrial);
        disp(count);
        disp(countIter);
   end
end 

cd(filePath)
save("avgMaxActivation.mat", "avgMaxActPerBirdArray");