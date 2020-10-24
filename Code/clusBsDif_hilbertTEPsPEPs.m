function [stat] = clusBsDif_hilbertTEPsPEPs(pathOut, pathIn, tepCond, freqToExamin, baseLine, dataName, latency)

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

% Load neighbour template
load([pathIn, 'neighbour_template.mat']);


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
    
    
% Data to fieldtrip format
for idx = 1:length(ID)
    allData.condA{idx} = [];
    allData.condA{idx}.label = eeglabChans;
    allData.condA{idx}.time = [-1:0.001:0.999];
    allData.condA{idx}.dimord = 'chan_time';
    allData.condA{idx}.avg = squeeze(freqDataA(:,:,idx));
    allData.condB{idx} = [];
    allData.condB{idx}.label = eeglabChans;
    allData.condB{idx}.time = [-1:0.001:0.999];
    allData.condB{idx}.dimord = 'chan_time';
    allData.condB{idx}.avg = squeeze(freqDataB(:,:,idx));
end


% % Load GrandAverage file
ft_defaults;
cfg = [];
cfg.channel = {'all'};
cfg.latency = latency;
cfg.method = 'montecarlo';
cfg.statistic = 'depsamplesT';
cfg.correctm = 'cluster';
cfg.clusteralpha = 0.05;
cfg.clusterstatistic = 'maxsum';
cfg.minnbchan = 2;
cfg.neighbours = neighbours;
cfg.tail = 0;
cfg.clustertail = 0;
cfg.alpha = 0.025;
cfg.numrandomization = 5000;
subj = length(ID);
design = zeros(1,2*subj);

for i = 1:subj
    design(1,i) = i;
end
for i = 1:subj
    design(1,subj+i) = i;
end
design(2,1:subj)        = 1;
design(2,subj+1:2*subj) = 2;
cfg.design = design;
cfg.uvar  = 1;
cfg.ivar  = 2;

% Clusterbased permutation test for the selected conditions
[stat] = ft_timelockstatistics(cfg,allData.condA{:}, allData.condB{:});
end

