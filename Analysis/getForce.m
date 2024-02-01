function [Force, Force_iso] = getForce(l_ce, v_ce, EMG, model_name, model_params, muscle_name, bird_name, bird_data)

if nargin > 6
    musvar = getMuscleParameters(bird_data, bird_name, muscle_name);
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

