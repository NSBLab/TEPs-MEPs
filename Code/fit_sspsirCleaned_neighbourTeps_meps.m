function [rho, randPerm_rho, pVal, pFDR] = fit_sspsirCleaned_neighbourTeps_meps(pathIn, tepCond, mepCond, Data, centChan, preStim, fl, constrainTo, nI, tFlip, winLength)
                                                                               
neighbour = load([pathIn,'neighbour_template.mat']);
Chans = neighbour.neighbours(strcmpi({neighbour.neighbours.label},centChan)).neighblabel;
Chans{end+1} = centChan;

% Load TEPs
load([pathIn,Data,'_',tepCond,'_TEPs.mat'],'suppr_data_SIR','ID','eeglabChans');

% Load MEPs
load([pathIn,'ppTMS_MEPs.mat'],'all_MEPs','mepCondition','ISIs');

mepType= strcmp(mepCondition,mepCond);
for jj = 1:length(Chans)
    j(jj) = find(strcmpi(eeglabChans,Chans{jj}));
end


if strcmp(fl,'non-flip')
    for idx = 1:length(ID)
        dtEx{idx} = double(squeeze(suppr_data_SIR(:,idx,:)));
        dataToExamin{idx} = mean(dtEx{idx}(j,:),1);
    end
elseif strcmp(fl,'flip')
     win1 = [1:preStim+tFlip];
     win2 = [length(win1):preStim+winLength-1];
    for idx = 1:length(ID)
        dtEx{idx} = double(squeeze(suppr_data_SIR(:,idx,:)));
        orig_meanTrials{idx} = mean(dtEx{idx}(j,:),1);
        a =[];
        b = [];
        a = -(orig_meanTrials{idx}(:,win1));
        b = (orig_meanTrials{idx}(:,win2));
        flip_meanTrials{idx} = [a b];
        ms(idx,:) = flip_meanTrials{idx};
        
    end
    flip_meanSubject = squeeze(mean(ms,1));
    dataToExamin = flip_meanTrials;
end

% Offset to zero-center EEG and EMG signals
EMG_offset = 100;
EEG_offset = 0;

% MEP side
MEP = all_MEPs{mepType}-EMG_offset;

% In case a peak is selected to be removed
if exist('pN','var')
    removeP = find(ISIs==pN);
    %   removeP = [1 2];
    ISIs(:,removeP)=[];
    MEP(:,removeP)=[];
end

peaks = preStim+ISIs;

% The selected TEP and MEP conditions
for idx = 1:length(ID)
    
    % Different TEP dataTypes
    TEP{idx} = dataToExamin{idx}-EEG_offset;
    TEP_Peaks{idx} = TEP{idx}(peaks);
    
    % Define your distance function inline measure of similarity between the two vectors as a function of the scaling parameter
    f_dist = @(eta)sum(abs(TEP_Peaks{idx}*eta-MEP(idx,:)));
    
    % Determine the best scaling parameter:
    if strcmp(constrainTo,'NegValues')
        eta_opt(idx) = fminsearchbnd (f_dist, 1,-inf,0);
    elseif strcmp(constrainTo,'PosValues')
        eta_opt(idx) = fminsearchbnd (f_dist, 1,0,inf);
    elseif strcmp(constrainTo,'Inf')
        eta_opt(idx) = fminsearchbnd (f_dist,1);
    end
    
    % Distance metric at best scaling parameter:
    rho(idx) = f_dist(eta_opt(idx));
    
    % Rescale teps
    rescaled_TEPs_peaks(idx,:) = TEP_Peaks{idx}*eta_opt(idx);
    rescaled_allTEPs(idx,:) = TEP{idx}*eta_opt(idx);
    
end

% Permuted distances for MEPs and TEPs acoss all conditions
for i = 1:nI
    rND_MEP = randperm(length(ISIs));
    rND_TEP = randperm(length(ISIs));
    
    for idx = 1:length(ID)
        
        % Random permutation of MEPs
        randPermMEPs{idx}(i,:) = MEP(idx,rND_MEP);
        
        % Random permutation of TEPs for the selected channel
        randPerm_TEP_Peaks{idx}(i,:) = TEP_Peaks{idx}(rND_TEP);
        
        % Calulates the minimum distance between randomly permuted TEPs and MEPs
        randPerm_f_dist = @(eta)sum(abs(randPerm_TEP_Peaks{idx}(i,:)*eta-randPermMEPs{idx}(i,:)));
        
        % Determine the best scaling parameter for each channel at each permutation iteration
        if strcmp(constrainTo,'NegValues')
            randPerm_eta_opt(idx,i) = fminsearchbnd (randPerm_f_dist, 1,-inf,0);
        elseif strcmp(constrainTo,'PosValues')
            randPerm_eta_opt(idx,i) = fminsearchbnd (randPerm_f_dist, 1,0,inf);
        elseif strcmp(constrainTo,'Inf')
            randPerm_eta_opt(idx,i) = fminsearchbnd (randPerm_f_dist,1);
        end
        
        % Distance metric at best scaling parameter for each channel at each permutation iteration
        randPerm_rho(idx,i) = randPerm_f_dist(randPerm_eta_opt(idx,i));
        
    end
end

sortPerm = sort(sum(randPerm_rho));
nPo = find(sortPerm <= sum(rho));
pVal =  max(nPo)/nI;
pFDR = mafdr(pVal,'BHFDR',true);

end