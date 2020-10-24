function [meanR Pval] = hilbertTEPsPEPs_temporalCorr(pathOut, pathIn, tepCond, freqToExamin, baseLine, dataName)

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
        freqDataA =  dataA.Alpha_bsLineCorrected;
        freqDataB = dataB.Alpha_bsLineCorrected;
    elseif strcmp(freqToExamin,'Beta')
        freqDataA =  dataA.Beta_bsLineCorrected;
        freqDataB = dataB.Beta_bsLineCorrected;
    elseif strcmp(freqToExamin,'Gamma')
        freqDataA =  dataA.Gamma_bsLineCorrected;
        freqDataB = dataB.Gamma_bsLineCorrected;
    elseif strcmp(freqToExamin,'BetaGamma')
        freqDataA =  dataA.BetaGamma_bsLineCorrected;
        freqDataB = dataB.BetaGamma_bsLineCorrected;
    end
elseif strcmp(baseLine,'uncorrected')
    if strcmp(freqToExamin,'Alpha')
        freqDataA =  dataA.Alpha;
        freqDataB = dataB.Alpha;
    elseif strcmp(freqToExamin,'Beta')
        freqDataA =  dataA.Beta;
        freqDataB = dataB.Beta;
    elseif strcmp(freqToExamin,'Gamma')
        freqDataA =  dataA.Gamma;
        freqDataB = dataB.Gamma;
    elseif strcmp(freqToExamin,'BetaGamma')
        freqDataA =  dataA.BetaGamma;
        freqDataB = dataB.BetaGamma;
    end
end

for idx= 1:length(ID)
    %correlation across conditions for Each electrode each subject
    for j = 1:length(eeglabChans)
        [Cor(j,idx), Pval(j,idx)] = corr(freqDataA(j,:,idx)' , freqDataB(j,:,idx)', 'Type','Spearman');
    end
end

%fisher's r to z transformation
z = .5.*log((1+Cor)./(1-Cor));

% Meansubjects from z and replace Inf with 1 and -inf with -1
meanz = mean(z,2) ;
infPos=[];
infNeg=[];
infPos= find(meanz==inf);
meanz(infPos)= 1;
infNeg= find(meanz==-inf);
meanz(infNeg) = -1;

% One sample permutaion test
[Pval] = mult_comp_perm_t1(z',50000,1);

% Transform meansubjects' z back to r for plotting
for j = 1:length(eeglabChans)
    meanR(j,1) = (exp(1)^(2.*meanz(j))-1)/(exp(1)^(2.*meanz(j))+1);
end
end

