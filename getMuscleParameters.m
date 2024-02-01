function musvar = getMuscleParameters(bird_data, bird_name, muscle_name)

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
