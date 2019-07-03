% EEGLAB, TESA and FieldTrip should be added to the path

% Each section has three parts : 1) Define options, 2) Run the function, and 3) Save the outputs

% 1) Define options : the chosen options here are specific to the current study taking the data from the Input folder and saving the results to the Output folder.
% The options with stars indicate that more than one choices are available for this study (choices are also provided).

% Input folder should contain:

% 1) ppTMS_MEPs.mat including:
% a) ppTMS-MEPs data (all_MEPs)
% b) Conditioned MEPs(all_condMEPs)
% c) Test (unconditioned) MEPs(all_testMEPs)
% d) ISIs for ppTMS (ISIs)
% e) Name of MEP conditions (mepCondition)
% f) Name of test MEPs conditions (testMepCondition)


% 2) TEP.mat and/or SSPSIR-Filtered-TEPs.mat including: 
% a) Name of TEP conditions (condition)
% b) EEGLAB channel file (eeglabChans)
% c) IDs of subjects (ID)
% d) TEP data for each subject at each condition (meanSubject)
% e) Mean of TEPs across subjects for each condition (meanTrials)

% For this study download the Input folder from https://doi.org/10.26180/5c0c8bf85eb24 and Unzip the folder
% ID of each experiment ('Biphasic' or 'Monophasic') is indicated in each file's name. What you need to do is:
% 1) Remove the _ID from the MEP file's name
% 2) Make a folder within Inputs and Outputs folders for each ID
% 3) Name it as the ID of the experiment ('Biphasic' or 'Monophasic' )
% 4) Save the MEP files to their corresponding folder
% 5) Move this script (preprocessingPipleline_MB.m) and AddPaths.m to the root path

% 2) Run the function: simply runs the function to do the nalysis for each section using the specified options.
% The functions can be used for other similar studies using the study-specific variables as "Options".

% 3) Save the outputs: saves the outputs of each function to the Output folder whitin the folder of the corresponding experimental condition ('Biphasic' or 'Monophasic')

%  Copyright (C) 2019, Mana Biabani <Mana.biabani@gmail.com> Monash University

clear; close all; clc;
AddPaths;

%%
% --------------------------------------------------------------------------------------------------
% ppTMS-MEPs : testing the size of changes in TS-induced MEPs caused by CS in ppTMS
% --------------------------------------------------------------------------------------------------

clear; close all; clc;

% Define the experimental condition (Name/ID of the experiment)
ExpCond = 'Monophasic'; % options : 'Biphasic','Monophasic' ***

% Where you get MEPs from
pathIn = ([pwd,'/Inputs/',ExpCond,'/']);

% Where you you save the outputs to
pathOut = ([pwd,'/Outputs/',ExpCond,'/']);

% % Load MEPs
% load([pathIn 'ppTMS_MEPs.mat']);

% Define conditioned MEPs
mepCond = 'left'; % Options for Biphasic = low, high; Options for Monophasic = right,left  ***

% Define unconditioned MEPs
testMepCond = 'left'; % Options for Biphasic = left; Options for Monophasic = right,left  ***
%---------------------------------------------------------------------------------------------------

% Run the function
[pVal, pFDR] = wilcoxonMEPs(pathIn, mepCond,testMepCond);

% Save the outputs
save([pathOut,'mepCond_',mepCond,'_testMepCond_',testMepCond,'_wilcoxonMEPs.mat']);

%%
% --------------------------------------------------------------------------------------------------
% Correlation between TEPs and MEPs amplitude
% --------------------------------------------------------------------------------------------------
clear; close all; clc;

% Experiment ID
ExpCond = 'Biphasic'; % options : 'Biphasic','Monophasic' ***

% Where you get MEPs and TEPs from
pathIn = ([pwd,'/Inputs/',ExpCond,'/']);

% Where you you save the outputs to
pathOut = ([pwd,'/Outputs/',ExpCond,'/']);

% Name of the .mat file containing TEPs
Data = 'TEPs';

% Which condition you want to get TEPs from
tepCond = 'high'; % Options for Biphasic = low, high; Options for Monophasic = right,left  ***

% Which condition you want to get MEPs from
mepCond ='high'; % Options for Biphasic = low, high; Options for Monophasic = right,left  ***

% Define the channel that you want to pick the neighbours of
centChan = 'C3';

% Define the length of preStim duration
preStim = 1000;
%---------------------------------------------------------------------------------------------------

% Run the function
[mepTepCor, mepTepPval, pFDR] = tepsMepsCorr(pathIn, Data, tepCond, mepCond, centChan, preStim);

% Save
save([pathOut,tepCond,'_',mepCond,'_',Data,'_tepsMepsCorr.mat']);
%%
% --------------------------------------------------------------------------------------------------
% Correlation between the amplitude of MEPs and slope of changes in TEPs at each peak
% --------------------------------------------------------------------------------------------------

clear; close all; clc;

% Experiment ID
ExpCond = 'Biphasic'; % options : 'Biphasic','Monophasic' ***

% Where you get MEPs and TEPs from
pathIn = ([pwd,'/Inputs/',ExpCond,'/']);

% Where you you save the outputs to
pathOut = ([pwd,'/Outputs/',ExpCond,'/']);

% Name of the .mat file containing TEPs
Data = 'TEPs';

% Which condition you want to get TEPs from
tepCond = 'high'; % Options for Biphasic = low, high; Options for Monophasic = right,left  ***

% Which condition you want to get MEPs from
mepCond ='high'; % Options for Biphasic = low, high; Options for Monophasic = right,left  ***

% Define the channel that you want to pick the neighbours of
centChan = 'C3';

% Define the length of preStim duration
preStim = 1000;

% Define the duration of TEPs preceding each peak to calculate slope from
slopeLength = 30;

% Define the length of early recordings to remove due to TMS artefact
DataToRemove = preStim+29;
%---------------------------------------------------------------------------------------------------

% Run
[pFDR, mepTepCor, mepTepPval] = tepSlope_MepAmplitude_corr(pathIn, Data, tepCond, mepCond, centChan, preStim, slopeLength, DataToRemove);

% Save 
save([pathOut,Data,'_',tepCond,'-teps_',mepCond,'-meps_Chan-',centChan,'_tepSlope_MepAmplitude_corr.mat']);
%%
% --------------------------------------------------------------------------------------------------
% Calculate similarity in shape between TEPs and MEPs (TEPs from neighbouring channels)
% --------------------------------------------------------------------------------------------------
clear; close all; clc;

% Experiment ID
ExpCond = 'Biphasic'; % Options : 'Biphasic','Monophasic' ***

% Where you get MEPs and TEPs from
pathIn = ([pwd,'/Inputs/',ExpCond,'/']);

% Where you you save the outputs to
pathOut = ([pwd,'/Outputs/',ExpCond,'/']);

% Choose TEPs : before or after suppressing PEPs with SSPSIR
Data = 'SSPSIR-Filtered-TEPs'; % Options : 'TEPs', 'SSPSIR-Filtered-TEPs'***

% Which condition you want to get TEPs from
tepCond = 'high'; % Options for Biphasic = low, high; Options for Monophasic = right,left  ***

% Which condition you want to get MEPs from
mepCond ='high'; % Options for Biphasic = low, high; Options for Monophasic = right,left  ***

% Define the channel that you want to pick the neighbours of
centChan = 'C3'; % Options : 'C3';'C4';'Cz';'FPz';'POz' ***

% Define the length of preStim duration
preStim = 100; % Options : 100 (for sspsir-filtered) , 1000 (for unfiltered TEPs) ***

% Define if you want to flip signal
fl = 'flip'; % Options : 'flip', 'non-flip' ***

% Determine if you want to confine the best scaling parameter to positive or negetive values (inf = infinity (no confines))
constrainTo = 'PosValues'; % Options : 'PosValues', 'NegValues', 'Inf' ***

% Number of iterations for permutation tests
nI = 1000;

% Length of TEP time window to examine
winLength = 250 ;

% The timepoint at which TEP signals are flippd (do not change if fl = nonflip)
tFlip = 60;
%---------------------------------------------------------------------------------------------------

% Run
[rho, randPerm_rho] = fit_neighbourTeps_meps (pathIn, Data, tepCond, mepCond, centChan, preStim, fl, constrainTo, nI, winLength, tFlip);

% Save
save([pathOut,Data,'_',tepCond,'-teps_',mepCond,'-meps_Chan-',centChan,'_',fl,'_fit_neighbourTeps_meps.mat']);
%%
% --------------------------------------------------------------------------------------------------
%  Clustering EEG electrodes based on the correlation between TEPs recorded by each pair of channels
% --------------------------------------------------------------------------------------------------
clear; close all; clc;

% Experiment ID
ExpCond = 'Biphasic'; % options : 'Biphasic','Monophasic' ***

% Where you get MEPs and TEPs from
pathIn = ([pwd,'/Inputs/',ExpCond,'/']);

% Where you you save the outputs to
pathOut = ([pwd,'/Outputs/',ExpCond,'/']);

% Choose TEPs : before or after suppressing PEPs with SSPSIR
Data = 'TEPs'; % Options : 'TEPs', 'SSPSIR-Filtered-TEPs'*** 

% Which condition you want to get TEPs from
tepCond = 'high'; % Options for Biphasic = low, high; Options for Monophasic = right,left  ***

% Define the length of preStim duration
preStim = 1000; % Options : 100 (for sspsir-filtered) , 1000 (for unfiltered TEPs) ***

% Length of TEP time window to examine
winLength = 250;

% Define a Percentage of Max Linkage as the Clustering Threshold
Pr = 0.6;
%---------------------------------------------------------------------------------------------------

% Run
[avrgSilhScore, sdSilhScore, chanGroupsInd] = corrBasedClustering(pathIn, Data, tepCond, preStim, winLength, Pr);

% Save
save([pathOut,Data,'_',tepCond,'_corrBasedClustering.mat']);
%%
% --------------------------------------------------------------------------------------------------
% Calculate similarity in shape between TEPs (recorded at each cluster) and MEPs
% --------------------------------------------------------------------------------------------------
clear; close all; clc;

% Experiment ID
ExpCond = 'Biphasic'; % options : 'Biphasic','Monophasic' ***

% Where you get MEPs and TEPs from
pathIn = ([pwd,'/Inputs/',ExpCond,'/']);

% Where you you save the outputs to
pathOut = ([pwd,'/Outputs/',ExpCond,'/']);

% Choose TEPs : before or after suppressing PEPs with SSPSIR
Data = 'TEPs'; % Options : 'TEPs', 'SSPSIR-Filtered-TEPs'*** 

% Which condition you want to get TEPs from
tepCond = 'high'; % Options for Biphasic = low, high; Options for Monophasic = right,left  ***

% Which condition you want to get MEPs from
mepCond ='high'; % Options for Biphasic = low, high; Options for Monophasic = right,left  ***

% Define the condition you want to get the clusters from
clusterCond = 'high';

% Define the length of preStim duration
preStim = 100; % Options : 100 (for sspsir-filtered) , 1000 (for unfiltered TEPs)

% Define if you want to flip signal
fl = 'flip'; % Options : 'flip', 'non-flip'  ***

% Determine if you want to confine the best scaling parameter to positive or negetive values (inf = infinity (no confines))
constrainTo = 'PosValues'; % Options : 'PosValues', 'NegValues', 'Inf'  ***

% Number of iterations for permutation tests
nI = 1000;

% Length of TEP time window to examine
winLength = 250 ;

% The timepoint at which TEP signals are flippd (do not change if fl = nonflip)
tFlip = 60;
%---------------------------------------------------------------------------------------------------

% Run
[rho, randPerm_rho, pVal, pFDR] = fit_clusterTeps_meps (pathIn, pathOut, Data, tepCond, mepCond, clusterCond, preStim, fl, constrainTo, nI, winLength, tFlip);

% Save
save([pathOut,Data,'_',tepCond,'-teps_',mepCond,'-meps_cluster_',clusterCond,'_',fl,'_fit_clusterTeps_meps.mat']);
%%
%---------------------------------------------------------------------------------------------------
% Comparison of flipping the polarity at different timepoints(ISIs)(TEPs from neighbouring channels)
%---------------------------------------------------------------------------------------------------

clear; close all; clc;

% Experiment ID
ExpCond = 'Biphasic'; % options : 'Biphasic','Monophasic' ***

% Where you get MEPs and TEPs from
pathIn = ([pwd,'/Inputs/',ExpCond,'/']);

% Where you you save the outputs to
pathOut = ([pwd,'/Outputs/',ExpCond,'/']);

% Choose TEPs
Data = 'TEPs';

% Which condition you want to get TEPs from
tepCond = 'high'; % Options for Biphasic = low, high; Options for Monophasic = right,left  ***

% Which condition you want to get MEPs from
mepCond ='high'; % Options for Biphasic = low, high; Options for Monophasic = right,left  ***

% Define the channel that you want to pick the neighbours of
centChan = 'C3'; % Options : 'C3';'C4';'Cz';'FPz';'POz' ***

% Define the length of preStim duration
preStim = 1000;

% Determine if you want to confine the best scaling parameter to positive or negetive values (inf = infinity (no confines))
constrainTo = 'PosValues'; % Options : 'PosValues', 'NegValues', 'Inf' ***

% Number of iterations for permutation tests
nI = 1000;
%---------------------------------------------------------------------------------------------------

% Run
[rho, randPerm_rho, pVal, pFDR] = timeSpecificFlip_fit_neighbourTeps_meps(pathIn, Data, tepCond, mepCond, centChan, preStim, constrainTo, nI);
                                                                       
% Save
save([pathOut,Data,'_',tepCond,'-teps_',mepCond,'-meps_Chan-',centChan,'_timeSpecificFlip_fit_neighbourTeps_meps.mat']);
%%
%---------------------------------------------------------------------------------------------------
% Comparison of flipping the polarity at different timepoints(ISIs)(TEPs from corr-based clusters)
%---------------------------------------------------------------------------------------------------
clear; close all; clc;

% Experiment ID
ExpCond = 'Biphasic'; % options : 'Biphasic','Monophasic' ***

% Where you get MEPs and TEPs from
pathIn = ([pwd,'/Inputs/',ExpCond,'/']);

% Where you you save the outputs to
pathOut = ([pwd,'/Outputs/',ExpCond,'/']);

% Name of the .mat file containing TEPs
Data = 'TEPs'; 

% Which condition you want to get TEPs from
tepCond = 'high'; % Options for Biphasic = low, high; Options for Monophasic = right,left  ***

% Which condition you want to get MEPs from
mepCond ='high'; % Options for Biphasic = low, high; Options for Monophasic = right,left  ***

% Define the condition you want to get the clusters from
clusterCond = 'high';


% Define the length of preStim duration
preStim = 1000;

% Determine if you want to confine the best scaling parameter to positive or negetive values (inf = infinity (no confines))
constrainTo = 'PosValues'; % Options : 'PosValues', 'NegValues', 'Inf'

% Number of iterations for permutation tests
nI = 1000;
%---------------------------------------------------------------------------------------------------

% Run
[rho, randPerm_rho, pVal, pFDR] = timeSpecificFlip_fit_clusterTeps_meps (pathIn, pathOut, Data, tepCond, mepCond, clusterCond, preStim,constrainTo, nI);

% Save
save([pathOut,Data,'_',tepCond,'-teps_',mepCond,'-meps_cluster_',clusterCond,'_timeSpecificFlip_fit_clusterTeps_meps.mat']);