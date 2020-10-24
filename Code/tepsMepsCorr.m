function [mepTepCor, mepTepPval, pFDR] = tepsMepsCorr(pathIn, Data, tepCond, mepCond, centChan, preStim)

    % load MEPs
    load([pathIn 'ppTMS_MEPs.mat']);

    % Load the neighbouring electrodes from the template taken from fieldtrip
    neighbour = load([pathIn,'neighbour_template.mat']);

    % Extract the neighbours of the defined centChan
    Chans = neighbour.neighbours(strcmpi({neighbour.neighbours.label},centChan)).neighblabel;

    % Add centChan
    Chans{end+1} = centChan;

    % Load MEPs
    load([pathIn 'ppTMS_MEPs.mat']);
    mCond = strcmp(mepCondition,mepCond);

    MEP = all_MEPs{mCond};
    peaks = preStim+ISIs;

    % Load TEPs at each ISI for the defined channel
    load([pathIn,Data,'.mat']);
    cond = strcmp(condition,tepCond);

    % Index of the selected channels
    j = zeros (1,length(Chans));
    for jj = 1:length(Chans)
        j(jj) = find(strcmp(eeglabChans,Chans{jj}));
    end

    % TEPs of the selected channels at each ISI
    TEP = zeros (length(Chans),length(ISIs));
    for idx = 1:length(ID)
        TEP(idx,:) = mean(meanTrials{cond}{idx}(j,peaks));
    end

    % TEPs-MEPs correlation at each ISI
    mepTepCor = zeros (1,length(ISIs));
    mepTepPval = zeros (1,length(ISIs));
    for p = 1:length(peaks)
        [mepTepCor(p), mepTepPval(p)] = corr(MEP(:,p) , TEP(:,p), 'Type','Spearman');
    end

    % FDR correction of pvalues across timepoints
    pFDR = mafdr(mepTepPval,'BHFDR',true);

end

