function [Force, Force_iso] = getForce(l_ce, v_ce, EMG, model_name, model_params, muscle_name, bird_name, bird_data)

if nargin > 6
    %% Load parameters from excel
    ind_row = findBird(bird_name);

    warning('The columns of the xlsx file are hard coded, please make sure that they match your version of MuscleMorphologyData.xlsx')

    if strcmpi(muscle_name, 'lg')
%         musvar.penn_ang = bird_data(ind_row,7);
        musvar.PCSA = bird_data(ind_row,8)/1000/1000; %converted to m
%         musvar.mmass = bird_data(ind_row,5)/1000;
        musvar.l_opt = bird_data(ind_row,6)/1000;
    elseif strcmpi(muscle_name, 'df')
%         musvar.penn_ang = bird_data(ind_row,13);
        musvar.PCSA = bird_data(ind_row,14)/1000/1000;
%         musvar.mmass = bird_data(ind_row,11)/1000;
        musvar.l_opt = bird_data(ind_row,12)/1000;
    else
        error('Incorrect muscle name')
    end

    sigma = 3.6e5; % Max muscle stress, from: S M Cox, K L Easton, M Cromie Lear, R L Marsh, S L Delp, J Rubenson, The Interaction of Compliance and Activation on the Force-Length Operating Range and Force Generating Capacity of Skeletal Muscle: A Computational Study using a Guinea Fowl Musculoskeletal Model, Integrative Organismal Biology, Volume 1, Issue 1, 2019, obz022, https://doi.org/10.1093/iob/obz022
    musvar.f_max = sigma*musvar.PCSA;
else
    musvar.f_max = 1;
    musvar.l_opt = 1;
end

if strcmpi(model_name, 'hill')
    if strcmpi(model_params(end-3:end), '.mat')

        load(model_params) %Change file name to change model parameters
        %assign optimized variables
        modelvar.v_max = optvar.Position(1);
        modelvar.Arel = optvar.Position(2);
        modelvar.gmax = optvar.Position(3);
        modelvar.W = optvar.Position(4);
        modelvar.kPEE = optvar.Position(5);
        modelvar.PEEslack = optvar.Position(6);
    else
        modelvar.v_max = 10;
        modelvar.PEEslack = 1.2;
        modelvar.gmax = 1.5;
        modelvar.kPEE = 1/musvar.l_opt^2;
        modelvar.Arel = 0.25;
        modelvar.W = 0.4;
    end

    %dependent parameters         
    modelvar.c3 = modelvar.v_max*modelvar.Arel*(modelvar.gmax - 1.)/(modelvar.Arel + 1);

    %% Calculate Hill forces
    if abs(std(l_ce)) < 1e-4 %means to check force-velocity relationship
        Force = findMuscleForce(modelvar, musvar, ones(size(v_ce)), v_ce, EMG, 'gvce');
    elseif abs(std(v_ce)) < 1e-4 %means to check force-length relationship including PEE
        Force = findMuscleForce(modelvar, musvar, l_ce, zeros(size(l_ce)), EMG, 'flcePEE');
    else
        Force = findMuscleForce(modelvar, musvar, l_ce, v_ce, EMG);
    end
elseif strcmpi(model_name, 'nn')
    load(model_params) %Change file name to change model parameters

    Mdl = NN.Mdl;
    data_mean = NN.data_mean;
    data_std = NN.data_std;

    if abs(std(l_ce)) < 1e-4 %means to check force-velocity relationship
        lce_norm = doDataNormalization(ones(size(v_ce)), data_mean.lce,data_std.lce);
        vce_norm = doDataNormalization(v_ce, data_mean.vce,data_std.vce);
        EMG_norm = doDataNormalization(EMG, data_mean.EMG,data_std.EMG);
    elseif abs(std(v_ce)) < 1e-4 %means to check force-length relationship
        lce_norm = doDataNormalization(l_ce, data_mean.lce,data_std.lce);
        vce_norm = doDataNormalization(zeros(size(l_ce)), data_mean.vce,data_std.vce);
        EMG_norm = doDataNormalization(EMG, data_mean.EMG,data_std.EMG);
    elseif abs(std(EMG)) < 1e-4 %means to check PEE relationship
        lce_norm = doDataNormalization(l_ce, data_mean.lce,data_std.lce);
        vce_norm = doDataNormalization(zeros(size(l_ce)), data_mean.vce,data_std.vce);
        EMG_norm = doDataNormalization(EMG, data_mean.EMG,data_std.EMG);
    else      
        lce_norm = doDataNormalization(l_ce, data_mean.lce,data_std.lce);
        vce_norm = doDataNormalization(v_ce, data_mean.vce,data_std.vce);
        EMG_norm = doDataNormalization(EMG, data_mean.EMG,data_std.EMG);
    end

    X = [lce_norm vce_norm EMG_norm];
    %% Calculate NN forces
    Force_norm = predict(Mdl, X);
    Force = doDataDenormalization(Force_norm, data_mean.force, data_std.force);
elseif strcmpi(model_name, 'nn_emg')
    load(model_params) %Change file name to change model parameters

    Mdl = NN.Mdl;
    data_mean = NN.data_mean;
    data_std = NN.data_std;

    EMG_norm = doDataNormalization(EMG, data_mean.EMG,data_std.EMG);
    
    X = EMG_norm;
    %% Calculate NN forces
    Force_norm = predict(Mdl, X);
    Force = doDataDenormalization(Force_norm, data_mean.force, data_std.force);
else
    error('incorrect model_name, should be Hill or NN')
end

if nargout > 1
    Force_iso = musvar.f_max;
end

