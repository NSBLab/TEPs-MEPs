function [rho, randPerm_rho, pVal, pFDR] = timeSpecificFlip_fit_clusterTeps_meps (pathIn, pathOut, Data, tepCond, mepCond, clusterCond, preStim,constrainTo, nI)

% Load TEPs
load([pathIn,Data,'.mat'],'condition','meanTrials','ID');

% Load MEPs
load([pathIn,'ppTMS_MEPs.mat'],'all_MEPs','mepCondition','ISIs');

% Load clusters
load([pathOut,Data,'_',clusterCond,'_corrBasedClustering.mat'],'chanGroupsInd');

% Offset to zero-center EEG and EMG signals
EMG_offset = 100;
EEG_offset = 0;

mepType = strcmp(mepCondition,mepCond);
MEP = all_MEPs{mepType}-EMG_offset;

tepType = strcmp(condition,tepCond);

for jj = 1:length(chanGroupsInd)
    for idx = 1:length(ID)
        ClusterData(idx,jj,:) = (mean(meanTrials{tepType}{idx}(chanGroupsInd{jj},:),1))';
    end
end
dataToExamin = ClusterData;

% In case a peak is selected to be removed
if exist('pN','var')
    removeP = find(ISIs==pN);
    %   removeP = [1 2];
    ISIs(:,removeP)=[];
    MEP(:,removeP)=[];
end

peaks = preStim+ISIs;

% All possible combinations of -1 and 1 in a vector with the same length as the peaks
A = ones(length(peaks));
B = tril(A);
nB = -B;
fB = nB==0;
nB(fB) = 1;
polarityOrder = [A(1,:);nB];

% Find the minimum distance between MEPs and TEPs with each polarity reversal window for each individual
for idx = 1:length(ID)
    
    % Pick TEPs at ISIs timepoints
    TEP{idx} = squeeze(dataToExamin(idx,:,:)-EEG_offset);
    TEP_Peaks{idx} = TEP{idx}(:,peaks);
    
    for po = 1:length(polarityOrder)
        newTep_Peaks=[];
        newTep_Peaks = TEP_Peaks{idx}.*polarityOrder(po,:);
        
        for clustNum = 1:length(chanGroupsInd)
            % Define your distance function inline measure of similarity between the two vectors as a function of the scaling parameter, eta (I used b in our meeting today):
            f_dist = @(eta)sum(abs(newTep_Peaks(clustNum,:)*eta-MEP(idx,:)));
            eta_opt = [];
            
            % Determine the best scaling parameter:
            if strcmp(constrainTo,'NegValues')
                eta_opt = fminsearchbnd (f_dist, 1,-inf,0);
            elseif strcmp(constrainTo,'PosValues')
                eta_opt = fminsearchbnd (f_dist, 1,0,inf);
            elseif strcmp(constrainTo,'Inf')
                eta_opt = fminsearchbnd (f_dist,1);
            end
            
            % Distance metric at best scaling parameter for each polarity reversal window:
            rho{po}(idx,clustNum) = f_dist(eta_opt);
            
            rhoOrder(idx,po,clustNum) = rho{po}(idx,clustNum);
            
            % The polarity order which gives the minimum distance for each subject at each cluster
            aa = [];
            aa = find(rhoOrder(idx,:,clustNum) == min(rhoOrder(idx,:,clustNum)));
            optOrderInd(idx,clustNum) = aa(1);
        end
    end
    
end

for clustNum = 1:length(chanGroupsInd)
    % The optimum order for each subject and each cluster
    optOrderEachIndiv(:,:,clustNum) = polarityOrder(optOrderInd(:,clustNum),:);
    
    % The most frequent optimum polarity for each peak across individuals
    chosenPolarityOrder(clustNum,:) = mode(optOrderEachIndiv(:,:,clustNum));
end

% Permuted distances for MEPs and each polarity order of TEPs  
for i = 1:nI
    rND_MEP = randperm(length(ISIs));
    rND_TEP = randperm(length(ISIs));
    
    for po = 1:length(polarityOrder)
        for idx = 1:length(ID)
            
            newTep_Peaks=[];
            newTep_Peaks = TEP_Peaks{idx}.*polarityOrder(po,:);
            
            % Random permutation of MEPs
            randPermMEPs = [];
            randPermMEPs = MEP(idx,rND_MEP);
            
            for clustNum = 1:length(chanGroupsInd)
                
                % Random permutation of TEPs for each cluster
                randPerm_TEP_Peaks = [];
                tpj = [];
                tpj = newTep_Peaks(clustNum,:);
                randPerm_TEP_Peaks = tpj(rND_TEP);
                
                % Calulates the minimum distance between randomly permuted TEPs (at each polarity order) and MEPs
                randPerm_f_dist = @(eta)sum(abs(randPerm_TEP_Peaks*eta-randPermMEPs));
                
                randPerm_eta_opt = [];
                % Determine the best scaling parameter for each channel at each permutation iteration and for each polarity order
                if strcmp(constrainTo,'NegValues')
                    randPerm_eta_opt = fminsearchbnd (randPerm_f_dist, 1,-inf,0);
                elseif strcmp(constrainTo,'PosValues')
                    randPerm_eta_opt = fminsearchbnd (randPerm_f_dist, 1,0,inf);
                elseif strcmp(constrainTo,'Inf')
                    randPerm_eta_opt = fminsearchbnd (randPerm_f_dist,1);
                end
                
                % Distance metric at best scaling parameter for each channel at each permutation iteration and for each polarity order
                randPerm_rho{po}(idx,clustNum,i) = randPerm_f_dist(randPerm_eta_opt);
                
            end
        end
    end
end

% Level of significance
for clustNum = 1:length(chanGroupsInd)
    for po = 1:length(polarityOrder)
        sortPerm = sort(sum(squeeze(randPerm_rho{po}(:,clustNum,:))));
        nPo = find(sortPerm <= sum(rho{po}(:,clustNum)));
        pVal(clustNum,po) =  max(nPo)/nI;
    end

    pFDR (clustNum,:) = mafdr(pVal(clustNum,:),'BHFDR',true);
end

end

