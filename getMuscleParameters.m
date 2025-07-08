function musvar = getMuscleParameters(bird_data, bird_name, muscle_name)

%% Load parameters from excel
ind_row = findBird(bird_name);

if strcmpi(muscle_name, 'lg')
%         musvar.penn_ang = bird_data(ind_row,7);
    musvar.PCSA = bird_data(ind_row,8)/1000/1000; %converted to m^2
%         musvar.mmass = bird_data(ind_row,5)/1000;
    musvar.l_opt = bird_data(ind_row,6)/1000;
elseif strcmpi(muscle_name, 'df')
%         musvar.penn_ang = bird_data(ind_row,13);
    musvar.PCSA = bird_data(ind_row,17)/1000/1000;
%         musvar.mmass = bird_data(ind_row,11)/1000;
    musvar.l_opt = bird_data(ind_row,15)/1000;
else
    error('Incorrect muscle name')
end

sigma = 2.879e5; % Max muscle stress, units=Pa or N/m^2, from Schwaner et al. "Linking in vivo muscle dynamics to force–length and force–velocity properties reveals that guinea fowl lateral gastrocnemius operates at shorter than optimal lengths"
musvar.f_max = sigma*musvar.PCSA;
