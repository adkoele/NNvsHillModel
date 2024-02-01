function obj = objective(x, modelvar, musvar, l_ce, v_ce, EMG, measuredForce_N)

modelvar.v_max = x(1);
modelvar.Arel = x(2);
modelvar.gmax = x(3);
modelvar.c3 = modelvar.v_max*modelvar.Arel*(modelvar.gmax - 1.)/(modelvar.Arel + 1);
modelvar.W = x(4);
modelvar.kPEE = x(5)/musvar.l_opt^2;
modelvar.PEEslack = x(6);

act = EMG;
    
Fsee = findMuscleForce(modelvar, musvar, l_ce, v_ce, act);

obj = sum((Fsee'*musvar.f_max-measuredForce_N).^2)/length(l_ce);

