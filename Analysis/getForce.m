function [Force, Force_iso] = getForce(folder_networksparams, l_ce, v_ce, EMG, model_name, model_params, muscle_name, bird_name, bird_data, nnMdlIdx)

if nargin > 7
    musvar = getMuscleParameters(bird_data, bird_name, muscle_name);
else
    musvar.f_max = 1;
    musvar.l_opt = 1;
end

if strcmpi(model_name, 'hill')
    if strcmpi(model_params(end-3:end), '.mat')
        load([folder_networksparams model_params])
        if contains(model_params, 'ModelwidthOptOnDF') && contains(model_params, 'r01')
            modelvar.v_max = 10.2103;
            modelvar.PEEslack = 2.7768; 
            modelvar.gmax = 0.9725;
            modelvar.kPEE = -1.8124;
            modelvar.Arel = 0.0278; 
            modelvar.W = optvar.Position(1); 
        elseif contains(model_params, 'ModelwidthOptOnDF') && contains(model_params, 'r12')
            modelvar.v_max = 4.5614;
            modelvar.PEEslack = 1.0002; 
            modelvar.gmax = 1.0427; 
            modelvar.kPEE = 0.3612;
            modelvar.Arel = 0.2263; 
            modelvar.W = optvar.Position(1); 
        else
        %assign optimized variables
            modelvar.v_max = optvar.Position(1);
            modelvar.Arel = optvar.Position(2);
            modelvar.gmax = optvar.Position(3);
            modelvar.W = optvar.Position(4);
            modelvar.kPEE = optvar.Position(5);
            modelvar.PEEslack = optvar.Position(6);
        end
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
elseif strcmpi(model_name, 'MillHill')
    if strcmpi(model_params(end-3:end), '.mat')
        load([folder_networksparams model_params])
        %assign optimized variables
        curveParams.v_max = optvar.Position(1);
        curveParams.eZero = optvar.Position(2);
        curveParams.eIsoY = optvar.Position(3);
        curveParams.kLowY = optvar.Position(4);
        curveParams.fpeeCurviness = optvar.Position(5);
        curveParams.fmaxE = optvar.Position(6);
        curveParams.dydxEY = optvar.Position(7);
        curveParams.dydxC = optvar.Position(8);
        curveParams.dydxNearEY = optvar.Position(9);
        curveParams.fvAtHalfVmax = optvar.Position(10);
        curveParams.gvceEccCurviness = optvar.Position(11);
        curveParams.lce0 = optvar.Position(12);
        curveParams.lce1Y = optvar.Position(13);
        curveParams.lce3Y = optvar.Position(14);
        curveParams.minActiveForceLengthValue = optvar.Position(15);
        curveParams.plateauSlopeY = optvar.Position(16);
        curveParams.flceCurviness = optvar.Position(17);   

        %dependent parameters
        eIsoMax = 1.2; % make sure that this value does not exist in the value range of eZero and is larger than that value range        
        curveParams.eIso = curveParams.eZero + curveParams.eIsoY * (eIsoMax - curveParams.eZero);
        curveParams.kIso = 2 / (curveParams.eIso - curveParams.eZero);
        curveParams.kLow = curveParams.kLowY * curveParams.kIso;
        curveParams.dydxE = curveParams.dydxEY * (curveParams.fmaxE-1);
        curveParams.dydxNearE = curveParams.dydxE + curveParams.dydxNearEY * (curveParams.fmaxE-1-curveParams.dydxE);
        curveParams.lce1 = curveParams.lce0 + curveParams.lce1Y;
        curveParams.lce2 = 1;
        curveParams.lce3 = curveParams.lce2 + curveParams.lce3Y;
        curveParams.plateauSlope = curveParams.plateauSlopeY * ((1-curveParams.minActiveForceLengthValue)/(curveParams.lce2-curveParams.lce1));

    else
        %If .mat file with optimized parameters hasn't been passed, then use
        %Millard's values (they correspond to human anatomy). 
        curveParams.v_max = 10;
        curveParams.eZero = 0;
        curveParams.eIsoY = 0.584;
        curveParams.kLowY = 0.07;
        curveParams.fpeeCurviness = 0.75;
        curveParams.fmaxE = 1.4;
        curveParams.dydxEY = 0.1;
        curveParams.dydxC = 0;
        curveParams.dydxNearEY = 0.17;
        curveParams.fvAtHalfVmax = 0.15;
        curveParams.gvceEccCurviness = 0.9;
        curveParams.lce0 = 0.45;
        curveParams.lce1Y = 0.28;
        curveParams.lce3Y = 0.8;
        curveParams.minActiveForceLengthValue = 0.1;
        curveParams.plateauSlopeY = 0.689;
        curveParams.flceCurviness = 0.9;

        %dependent params. here their values match those used by Millard in
        %his original code
        %dependent parameters
        eIsoMax = 1.2; % make sure that this value does not exist in the value range of eZero and is larger than that value range        
        curveParams.eIso = 0.7;
        curveParams.kIso = 2 / (curveParams.eIso - curveParams.eZero);
        curveParams.kLow = 0.2;
        curveParams.dydxE = 0.1;
        curveParams.dydxNearE = 0.15;
        curveParams.lce1 = 0.73;
        curveParams.lce2 = 1;
        curveParams.lce3 = 1.8123;
        curveParams.plateauSlope = 0.8616;

    end

    
    %% Calculate Millard Hill forces
    if abs(std(l_ce)) < 1e-4 %Calculate force to check force-velocity relationship
        
        optimizationFlag = false;
        v_ce = v_ce/10; 
        
        Force = findBezierCurvesAndMillardMuscleForce(curveParams, musvar, ones(size(v_ce)), v_ce, EMG, optimizationFlag, 'gvce');
    elseif abs(std(v_ce)) < 1e-4 %Calculate force to check force-length relationship including PEE
       
        optimizationFlag = false;
        Force = findBezierCurvesAndMillardMuscleForce(curveParams, musvar, l_ce, zeros(size(l_ce)), EMG, optimizationFlag, 'flcePEE');
    else
        v_ce = v_ce/curveParams.v_max;
        optimizationFlag = false;

        Force = findBezierCurvesAndMillardMuscleForce(curveParams, musvar, l_ce, v_ce, EMG, optimizationFlag);
    end

elseif strcmpi(model_name, 'nn')
    load([folder_networksparams model_params])

    if nargin<10 %means that no model index has been passed, so we will use the best NN
        Mdl = NN.Mdl;
    else    
        Mdl = res(nnMdlIdx).Mdl;
    end
    data_mean = NN.data_mean;
    data_std = NN.data_std;

    if abs(std(l_ce)) < 1e-4 %means to check force-velocity relationship
        lce_norm = doDataNormalization(ones(size(v_ce)), data_mean.lce,data_std.lce);
        vce_norm = doDataNormalization(v_ce, data_mean.vce,data_std.vce);
        EMG_norm = EMG; 
    elseif abs(std(v_ce)) < 1e-4 %means to check force-length relationship
        lce_norm = doDataNormalization(l_ce, data_mean.lce,data_std.lce);
        vce_norm = doDataNormalization(zeros(size(l_ce)), data_mean.vce,data_std.vce);
        EMG_norm = EMG; 
    elseif abs(std(EMG)) < 1e-4 %means to check PEE relationship
        lce_norm = doDataNormalization(l_ce, data_mean.lce,data_std.lce);
        vce_norm = doDataNormalization(zeros(size(l_ce)), data_mean.vce,data_std.vce);
        EMG_norm = EMG; 
    else      
        lce_norm = doDataNormalization(l_ce, data_mean.lce,data_std.lce);
        vce_norm = doDataNormalization(v_ce, data_mean.vce,data_std.vce);
        EMG_norm = EMG; 
    end

    X = [lce_norm vce_norm EMG_norm];
    %% Calculate NN forces
    Force_norm = predict(Mdl, X);
    Force = doDataDenormalization(Force_norm, data_mean.force, data_std.force);
elseif strcmpi(model_name, 'nn_emg')
    load([folder_networksparams model_params])

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

