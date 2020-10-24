function [rescaled_allTEPs, rho, randPerm_rho, pVal] = fit_hilbTepsMeps(pathOut, pathIn, mepCond, tepCond, freqToExamin, baseLine, dataName, centChan, constrainTo, preStim, nI)

% Name of the .mat file containing TEPs
load([pathOut,tepCond,'_hilbertEvokedAmp.mat'],dataName);
data = EvokedAmp;
clearvars  EvokedAmp

% Name of the .mat file containing PEPs
load([pathOut,'control_hilbertEvokedAmp.mat'],dataName);
dataB = EvokedAmp;
clearvars  EvokedAmp

% Load participant IDs
load([pathIn 'TEPs.mat'],'ID');

% Load eeg channels names (1 by channel cell)
load('eeglabChans.mat');

% Extract the neighbours of the defined centChan
load([pathIn 'neighbour_template.mat']);
Chan = neighbours(strcmpi({neighbours.label},centChan)).neighblabel;
Chan{end+1} = centChan;

% Index of the selected channels (to measure power from their average)
j = find(ismember(eeglabChans,Chan));

% Load MEPs
load([pathIn 'ppTMS_MEPs.mat']);
mepType= find(strcmp(mepCondition,mepCond));

% Offset to zero-center EEG and EMG signals
EMG_offset = 100;
EEG_offset = 0;

% MEP side
MEP = all_MEPs{mepType}-EMG_offset;

% Remove peak>10ms removed
if ismember(10, ISIs)
    removeP = find(ISIs==10);
    ISIs(:,removeP)=[];
    MEP(:,removeP)=[];
end

peaks = repmat((preStim+ISIs)',1, length(ID));

if strcmp(baseLine,'corrected')
    if strcmp(freqToExamin,'Alpha')
        freqData =  data.Alpha_bsLineCorrected;
    elseif strcmp(freqToExamin,'Beta')
        freqData =  data.Beta_bsLineCorrected;
    elseif strcmp(freqToExamin,'Gamma')
        freqData =  data.Gamma_bsLineCorrected;
    elseif strcmp(freqToExamin,'BetaGamma')
        freqData =  data.BetaGamma_bsLineCorrected;
    end
elseif strcmp(baseLine,'uncorrected')
    if strcmp(freqToExamin,'Alpha')
        freqData =  data.Alpha;
    elseif strcmp(freqToExamin,'Beta')
        freqData =  data.Beta;
    elseif strcmp(freqToExamin,'Gamma')
        freqData =  data.Gamma;
    elseif strcmp(freqToExamin,'BetaGamma')
        freqData =  data.BetaGamma;
    end
end

dataToExamin = squeeze(mean(freqData(j,:,:),1));

% The selected TEP and MEP conditions
for idx = 1:length(ID)
    
    % Different TEP dataTypes
    TEP{idx} = dataToExamin(:,idx)-EEG_offset;
    TEP_Peaks{idx} = TEP{idx}(peaks(:,idx));
    
    % Define your distance function inline ? measure of similarity between the two vectors as a function of the scaling parameter, eta (I used b in our meeting today):
    f_dist = @(eta)sum(abs(TEP_Peaks{idx}'*eta-MEP(idx,:)));
    
    % Determine the best scaling parameter:
    if strcmp(constrainTo,'NegValues')
        eta_opt(idx) = fminsearchbnd (f_dist, 1,-inf,0);
    elseif strcmp(constrainTo,'PosValues')
        eta_opt(idx) = fminsearchbnd(f_dist, 1,0,inf);
    end
    
    % Distance metric at best scaling parameter:
    rho(idx) = f_dist(eta_opt(idx));
    
    % Rescale teps
    rescaled_TEPs_peaks(idx,:) = TEP_Peaks{idx}*eta_opt(idx);
    rescaled_allTEPs(idx,:) = TEP{idx}*eta_opt(idx);
    
end

% Permuted distances for MEPs and TEPs acoss all conditions
for i = 1:nI
    rND_MEP = randperm(length(ISIs));
    rND_TEP = randperm(length(ISIs));
    
    for idx = 1:length(ID)
        
        % Random permutation of MEPs
        randPermMEPs{idx}(i,:) = MEP(idx,rND_MEP);
        
        % Random permutation of TEPs for the selected channel
        randPerm_TEP_Peaks{idx}(i,:) = TEP_Peaks{idx}(rND_TEP);
        
        % Calulates the minimum distance between randomly permuted TEPs and MEPs
        randPerm_f_dist = @(eta)sum(abs(randPerm_TEP_Peaks{idx}(i,:)*eta-randPermMEPs{idx}(i,:)));
        
        % Determine the best scaling parameter for each channel at each permutation iteration
        if strcmp(constrainTo,'NegValues')
            randPerm_eta_opt(idx,i) = fminsearchbnd (randPerm_f_dist, 1,-inf,0);
        elseif strcmp(constrainTo,'PosValues')
            randPerm_eta_opt(idx,i) = fminsearchbnd (randPerm_f_dist, 1,0,inf);
            
        end
        
        % Distance metric at best scaling parameter for each channel at each permutation iteration
        randPerm_rho(idx,i) = randPerm_f_dist(randPerm_eta_opt(idx,i));
        
    end
end
sortPerm = sort(sum(randPerm_rho));
distr = find(sortPerm <= sum(rho));
pVal =  max(distr)/nI;

end