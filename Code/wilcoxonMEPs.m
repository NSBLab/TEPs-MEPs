function [pVal, pFDR] = wilcoxonMEPs(pathIn, mepCond,testMepCond)

% Load MEPs
load([pathIn 'ppTMS_MEPs.mat']);

condMEPs = all_condMEPs{strcmpi(mepCondition,mepCond)};

uncondMEPs = all_testMEPs{strcmpi(testMepCondition,testMepCond)};

pVal = zeros(length(ISIs),1);
norm = zeros(length(ISIs),1);

    for I = 1:length(ISIs)
        [pVal(I)] = signrank(condMEPs(:,I),uncondMEPs);
        norm(I) = kstest(condMEPs(:,I));
    end

pFDR = mafdr(pVal,'BHFDR',true);
end