function [stat] = clustBsDif_medianSplitTeps(pathOut,pathIn,tepCond,latency)

% Load median split TEPs in filedtrip format
load([pathOut,tepCond,'_medianSplitTeps.mat']);

% Load neighbour template
load([pathIn,'neighbour_template.mat']);

% Load participant IDs
load([pathIn 'TEPs.mat'],'ID');

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
cfg.neighbours = neighbours;  % same as defined for the between-trials experiment
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
design(2,1:subj) = 1;
design(2,subj+1:2*subj) = 2;
cfg.design = design;
cfg.uvar  = 1;
cfg.ivar  = 2;

% clusterbased permutation test for the selected conditions
[stat] = ft_timelockstatistics(cfg,medianSplitTeps.upperMed{:}, medianSplitTeps.lowerMed{:});

end
