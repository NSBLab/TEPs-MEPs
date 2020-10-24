function [Cor ConfidIntval] = hilbertTEPsPEPs_spatialCorr(pathOut, pathIn, tepCond, freqToExamin, baseLine, dataName, win)

% Name of the .mat file containing TEPs
load([pathOut,tepCond,'_hilbertEvokedAmp.mat'],dataName);
dataA = EvokedAmp;
clearvars  EvokedAmp

% Name of the .mat file containing PEPs
load([pathOut,'control_hilbertEvokedAmp.mat'],dataName);
dataB = EvokedAmp;
clearvars  EvokedAmp

% Load participant IDs
load([pathIn 'TEPs.mat'],'ID');

% Load eeg channels names (1 by channel cell)
load('eeglabChans.mat');

if strcmp(baseLine,'corrected')
    if strcmp(freqToExamin,'Alpha')
        freqDataA =  dataA.Alpha_bsLineCorrected(:,win,:);
        freqDataB = dataB.Alpha_bsLineCorrected(:,win,:);
    elseif strcmp(freqToExamin,'Beta')
        freqDataA =  dataA.Beta_bsLineCorrected(:,win,:);
        freqDataB = dataB.Beta_bsLineCorrected(:,win,:);
    elseif strcmp(freqToExamin,'Gamma')
        freqDataA =  dataA.Gamma_bsLineCorrected(:,win,:);
        freqDataB = dataB.Gamma_bsLineCorrected(:,win,:);
    elseif strcmp(freqToExamin,'BetaGamma')
        freqDataA =  dataA.BetaGamma_bsLineCorrected(:,win,:);
        freqDataB = dataB.BetaGamma_bsLineCorrected(:,win,:);
    end
elseif strcmp(baseLine,'uncorrected')
    if strcmp(freqToExamin,'Alpha')
        freqDataA =  dataA.Alpha(:,win,:);
        freqDataB = dataB.Alpha(:,win,:);
    elseif strcmp(freqToExamin,'Beta')
        freqDataA =  dataA.Beta(:,win,:);
        freqDataB = dataB.Beta(:,win,:);
    elseif strcmp(freqToExamin,'Gamma')
        freqDataA =  dataA.Gamma(:,win,:);
        freqDataB = dataB.Gamma(:,win,:);
    elseif strcmp(freqToExamin,'BetaGamma')
        freqDataA =  dataA.BetaGamma(:,win,:);
        freqDataB = dataB.BetaGamma(:,win,:);
    end
end

for idx= 1:length(ID)
    % correlation of topoplots between conditions at each point of time and for each subject
    for t = 1:length(win)
        [Cor(idx,t), Pval(idx,t)] = corr(freqDataA(:,t,idx) , freqDataB(:,t,idx), 'Type','Spearman');
    end
end

% fisher transformation of r values and calculate Confidence intervals of correlations across subjects
for t = 1:length(win)
    for idx = 1:length(ID)
        r =[];
        r(idx,t)= Cor(idx,t);
        %fisher's r to z transformation
        z(idx,t)=.5.*log((1+r(idx,t))./(1-r(idx,t)));
    end
    % CI of z scores
    CI(:,t)= confidence_intervals(z(:,t),95);
    % average of z scores
    avZ(t) = mean(z(:,t),1);
    % fisher's z to r tranformation
    rFromZ(t) = (exp(1)^(2.*avZ(t))-1)/(exp(1)^(2.* avZ(t))+1);
    % fisher's z to r tranformation for CI
    rCIrFromZ(1,t) = (exp(1)^(2.*CI(1,t))-1)/(exp(1)^(2.* CI(1,t))+1);
    rCIrFromZ(2,t) =(exp(1)^(2.*CI(2,t))-1)/(exp(1)^(2.* CI(2,t))+1);
end
Cor = rFromZ;
ConfidIntval = rCIrFromZ;
end

