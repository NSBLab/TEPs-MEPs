function [Cor, Pval, pFDR] = mepsHilbertTEPsCorr(pathOut, pathIn, mepCond, tepCond, freqToExamin, baseLine, dataName, centChan, preStim)

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

mCond = find(strcmp(mepCondition,mepCond));
MEP = all_MEPs{mCond};

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

tPoint = ISIs+preStim;

% Hilbert TEPs-MEPs correlation at each ISI
for p = 1:length(ISIs)
    [Cor(p), Pval(p)] = corr(MEP(:,p), squeeze(mean(freqData(j,tPoint(p),:),1)), 'Type','Spearman');
end

% FDR correction of pvalues across timepoints
pFDR = mafdr(Pval,'BHFDR',true);

end

