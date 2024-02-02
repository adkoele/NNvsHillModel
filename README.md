# NNvsHillModel
Code used for: Neural Networks Estimate Muscle Force in Dynamic Conditions Better than Hill-type Muscle Models

This repository contains four folders and some general files that are contained in the main folder. Note that the function RMSE is available in MATLAB from version R2022b onwards. You can download a similar function here: https://de.mathworks.com/matlabcentral/fileexchange/21383-rmse

Folder Analysis: contains the code that was used to analyse our neural networks and Hill-type muscle models. It contains the following functions:
- Analysis_singletrialNNvsHill: this function was used to compare the "smaller" neural networks agains the Hill-type muscle models that were optimized on the same data
- Analysis_multtialNN: this function was used to evaluate the neural network that was generated with a large dataset.
- getForce: this function is called by both analysis functions and calculates the muscle force using one of the neural networks or Hill-type muscle models.
- rsquare: this function calculates the coefficient of determination (R^2).
  
**IMPORTANT**: please adjust the path in lines 7, 16, and 17 in both Analysis_singletrialNNvsHill and Analysis_multtialNN to reflect your system's path.

Folder Hill_optimization: contains the code that was used to optimize the Hill-type muscle models.
- guinneaFowls_optHillModel: main function that runs the optimization (note that line 6 and 14 need to be adapted to your own folder structure).
- cmaes: the function used for the optimization using the CMA-ES algorithms.
- objective: the objective that was minimized in the optimization (comparison of Hill force to measured force).

Folder NN_training_singletrial: contains the code that was used to train the neural networks that were compared to the Hill-type muscle models.
- NN_training_singletrial: contains the code that was used to train the neural networks.
- NN_training_singletrial_EMGonly: contains a similar code that was used to a neural network that only used activation as input, and not muscle length and velocity.

Folder NN_training_multipletrials: contains the code that was used to train a neural network on a large dataset of guinea fowl data.
- NN_training_multtials: main file used for training (note that line 11 and line 18 need to be adapted to your own folder structure).
- createDataset.m: this function should be run first and creates the dataset that is loaded in line 21 of NN_training_multtrials.

Then, the following functions are saved in the main folder, since these are used both for training/optimization and analysis.
- doDataDenormalization: converts the normalized output of the neural network (mean 0 and standard deviation 1) to a force normalized to isometric force.
- doDataNormalization: converts the inputs (activation, length and velocity normalized to optimal fibre length) to normalized inputs (mean 0 and standard deviation 1) as they are used in the neural network.
- findBird: outputs the correct row in the excel file given the bird name.
- findMuscleForce: finds the force in the muscle using the Hill-type muscle model.
- getMuscleParameters: loads the parameters of the current muscle and bird from the (already loaded) MuscleMorphologyData.xlsx file
- loadDataFile: loads the guinea fowl data and performs the filtering and processing on the shortening velocity and electromyography signal.
