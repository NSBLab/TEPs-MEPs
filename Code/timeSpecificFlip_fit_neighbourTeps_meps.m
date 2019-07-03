function [rho, randPerm_rho, pVal, pFDR] = timeSpecificFlip_fit_neighbourTeps_meps (pathIn, Data, tepCond, mepCond, centChan, preStim, constrainTo, nI)

neighbour = load([pathIn,'neighbour_template.mat']);
Chans = neighbour.neighbours(strcmpi({neighbour.neighbours.label},centChan)).neighblabel;
Chans{end+1} = centChan;

% Load TEPs
load([pathIn,Data,'.mat'],'condition','meanTrials','ID','eeglabChans');

% Load MEPs
load([pathIn,'ppTMS_MEPs.mat'],'all_MEPs','mepCondition','ISIs');

for jj = 1:length(Chans)
    j(jj) = find(strcmpi(eeglabChans,Chans{jj}));
end

% Offset to zero-center EEG and EMG signals
EMG_offset = 100;
EEG_offset = 0;

mepType = strcmp(mepCondition,mepCond);
MEP = all_MEPs{mepType}-EMG_offset;

tepType = strcmp(condition,tepCond);

for idx = 1:length(ID)
    NeighbouringData(idx,:) = (mean(meanTrials{tepType}{idx}(j,:),1));
end

dataToExamin = NeighbouringData;

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
    TEP{idx} = squeeze(dataToExamin(idx,:)-EEG_offset);
    TEP_Peaks{idx} = TEP{idx}(:,peaks);
    
    for po = 1:length(polarityOrder)
        newTep_Peaks=[];
        newTep_Peaks = TEP_Peaks{idx}.*polarityOrder(po,:);
        
        for clustNum = 1
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
            
        end
    end
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
            
            for clustNum = 1
                
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
for po = 1:length(polarityOrder)
    sortPerm = sort(sum(squeeze(randPerm_rho{po})));
    nPo = find(sortPerm <= sum(rho{po}));
    pVal(po) =  max(nPo)/nI;
end
pFDR = mafdr(pVal,'BHFDR',true);

end