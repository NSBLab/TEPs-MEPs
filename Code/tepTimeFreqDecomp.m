function [dataOut] = tepTimeFreqDecomp(pathIn, Data, tepCond, baseLine, freqInt)
ft_defaults;
load([pathIn Data,'.mat']);

for idx = 1:length(ID)
    
    % Load data in FieldTrip format
    data     = [];
    cfg            = [];
    cond = strcmp(condition,tepCond);
    ftdata = data{cond}{idx};
    data = ft_timelockanalysis(cfg,ftdata);
    
    % Time-frequency analysis - wavelet
    dataWav        = [];
    cfg            = [];
    cfg.channel    = {'all'};
    cfg.output     = 'pow';
    cfg.method     = 'wavelet';
    cfg.width      = 3.5;
    cfg.foi = freqInt;
    cfg.toi = -1:0.0001:0.999;
    cfg.keeptrials = 'no';
    dataWav = ft_freqanalysis(cfg,data);
    
    % Normalised to baseline
    cfg =[];
    cfg.baseline = baseLine;
    cfg.baselinetype = 'db';
    dataWavBc = ft_freqbaseline(cfg,dataWav);
    
    % Data name to save
    newId = ['Subj',ID{idx}];
    dataOut.(newId) = dataWavBc;
end
end
