function [rho, randPerm_rho, pVal, pFDR] = fit_clusterTeps_meps (pathIn, pathOut, Data, tepCond, mepCond, clusterCond, preStim, fl, constrainTo, nI, winLength, tFlip)

% Offset to zero-center EEG and EMG signals
EMG_offset = 100;
EEG_offset = 0;

% Load MEPs
load([pathIn,'ppTMS_MEPs.mat'],'all_MEPs','mepCondition','ISIs');
mepType = strcmp(mepCondition,mepCond);
MEP = all_MEPs{mepType} - EMG_offset;

% Load TEPs
load([pathIn, Data,'.mat'],'condition','ID','meanTrials','eeglabChans');
tepType = strcmp(condition,tepCond);

% Electrodes in each cluster
load([pathOut,Data,'_',tepCond,'_corrBasedClustering.mat'],'chanGroupsInd');

dataToExamin = [];
a =[];
b =[];

for jj = 1:length(chanGroupsInd)
    for idx = 1:length(ID)
        ClusterData(idx,jj,:) = (mean(meanTrials{tepType}{idx}(chanGroupsInd{jj},:),1))';% exclude control
    end
end

if strcmp(fl,'non-flip')
    dataToExamin = ClusterData;
elseif strcmp(fl,'flip')
    win1 = [1:preStim+tFlip];
    win2 = [length(win1):preStim+winLength-1];
    a = -(ClusterData(:,:,win1));
    b = ClusterData(:,:,win2);
    dataToExamin = cat(3,a,b);
end

% In case a peak is selected to be removed
if exist('pN','var')
    removeP = find(ISIs == pN);
    %   removeP = [1 2];
    ISIs(:,removeP) = [];
    MEP(:,removeP) = [];
end
peaks = preStim + ISIs;


% The selected TEP and MEP conditions
for idx = 1:length(ID)
    
    % Different TEP dataTypes
    TEP{idx} = squeeze(dataToExamin(idx,:,:)-EEG_offset);
    TEP_Peaks{idx} = TEP{idx}(:,peaks);
    
    for j = 1:length(chanGroupsInd)
        % Define your distance function inline ? measure of similarity between the two vectors as a function of the scaling parameter, eta (I used b in our meeting today):
        f_dist = @(eta)sum(abs(TEP_Peaks{idx}(j,:)*eta-MEP(idx,:)));
        
        % Determine the best scaling parameter:
        if strcmp(constrainTo,'NegValues')
            eta_opt(idx,j) = fminsearchbnd (f_dist, 1,-inf,0);
        elseif strcmp(constrainTo,'PosValues')
            eta_opt(idx,j) = fminsearchbnd (f_dist, 1,0,inf);
        elseif strcmp(constrainTo,'Inf')
            eta_opt(idx,j) = fminsearchbnd (f_dist,1);
        end
        
        % Distance metric at best scaling parameter:
        rho(idx,j) = f_dist(eta_opt(idx,j));
        
        % Rescale teps
        rescaled_TEPs_peaks(idx,j,:) = TEP_Peaks{idx}(j,:)*eta_opt(idx,j);
        rescaled_allTEPs(idx,j,:) = TEP{idx}(j,:)*eta_opt(idx,j);
    end
    
end

% Permuted distances for MEPs and TEPs acoss all conditions
for i = 1:nI
    rND_MEP = randperm(length(ISIs));
    rND_TEP = randperm(length(ISIs));
    
    for idx = 1:length(ID)
        
        % Random permutation of MEPs
        randPermMEPs = [];
        randPermMEPs = MEP(idx,rND_MEP);
        
        for j = 1:length(chanGroupsInd)
            
            % Random permutation of TEPs for each channel
            randPerm_TEP_Peaks = [];
            tpj = [];
            tpj = TEP_Peaks{idx}(j,:);
            randPerm_TEP_Peaks = tpj(rND_TEP);
            
            % Calulates the minimum distance between randomly permuted TEPs and MEPs
            randPerm_f_dist = @(eta)sum(abs(randPerm_TEP_Peaks*eta-randPermMEPs));
            
            % Determine the best scaling parameter for each channel at each permutation iteration
            if strcmp(constrainTo,'NegValues')
                randPerm_eta_opt(idx,j,i) = fminsearchbnd (randPerm_f_dist, 1,-inf,0);
            elseif strcmp(constrainTo,'PosValues')
                randPerm_eta_opt(idx,j,i) = fminsearchbnd (randPerm_f_dist, 1,0,inf);
            elseif strcmp(constrainTo,'Inf')
                randPerm_eta_opt(idx,j,i) = fminsearchbnd (randPerm_f_dist,1);
            end
            
            % Distance metric at best scaling parameter for each channel at each permutation iteration
            randPerm_rho(idx,j,i) = randPerm_f_dist(randPerm_eta_opt(idx,j,i));
            
        end
    end
end

for k = 1:length(chanGroupsInd)
    sortPerm = sort(sum(squeeze(randPerm_rho(:,k,:))));
    nPo = find(sortPerm <= sum(rho(:,k)));
    pVal(k) =  max(nPo)/nI;
end
pFDR = mafdr(pVal,'BHFDR',true);

end

