function [pFDR, mepTepCor, mepTepPval] = tepSlope_MepAmplitude_corr (pathIn, Data, tepCond, mepCond, centChan, preStim, slopeLength, DataToRemove)

% Load the neighbouring electrodes from the template taken from fieldtrip
neighbour = load([pathIn,'neighbour_template.mat']);

% Extract the neighbours of the defined centChan
Chans = neighbour.neighbours(strcmpi({neighbour.neighbours.label},centChan)).neighblabel;

% Add centChan
Chans{end+1} = centChan;

% Load MEPs
load([pathIn 'ppTMS_MEPs.mat']);
mCond = strcmp(mepCondition,mepCond);
MEP = all_MEPs{mCond}-100;
peaks = preStim+ISIs;

% Load TEPs at each ISI for the defined channel
load([pathIn, Data,'.mat']);
cond = strcmp(condition,tepCond);

% Index of the selected channels
for jj = 1:length(Chans)
    j(jj) = find(strcmp(eeglabChans,Chans{jj}));
end

if exist ('DataToRemove', 'var')
    sP = find (peaks>DataToRemove);
    slopePeaks = peaks(sP);
    slopeMEPs = MEP(:,sP);
else
    slopePeaks = peaks;
    slopeMEPs = MEP; 
end

for p = 1:length(slopePeaks)
    slopeTime(p,:) = slopePeaks(p)-slopeLength:slopePeaks(p);
end

% TEPs of the selected channels at each duration for slope analysis
for p = 1:length(slopePeaks) 
    for idx = 1:length(ID)
        TEP{p}(idx,:) = mean(meanTrials{cond}{idx}(j,slopeTime(p,:)));
        derv(idx,p) = mean(diff(TEP{p}(idx,:)));
    end
end

% TEPs-MEPs correlation at each ISI
for p = 1:length(slopePeaks)
    [mepTepCor(p), mepTepPval(p)] = corr(slopeMEPs(:,p) , derv(:,p), 'Type','Spearman');
end

% FDR correction of pvalues across timepoints
pFDR = mafdr(mepTepPval,'BHFDR',true);
end