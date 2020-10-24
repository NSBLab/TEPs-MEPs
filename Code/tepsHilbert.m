function [EvokedAmp] = tepsHilbert(pathIn,ID,tepCond,baseLine,eeglabChans)

for idx = 1:length(ID)
    myEEG = [];
    myEEG = pop_loadset([pathIn, ID{idx,1},'_',tepCond, '_laplacian_FINAL.set']);
    % --------------------- filters on meantrials for the amplitude of evoked potentials ---------------
    % Filter and Hilbert transformation for each electrode
    for j = 1:length(eeglabChans)
        
        % Bandpass filter data
        EEG_Alpha_meanTrials = [];
        EEG_Beta_meanTrials = [];
        EEG_Gamma_meanTrials = [];
        EEG_BetaGamma_meanTrials = [];
        
        % Get the average of TEPs across trials for one channel
        newEEG = [];
        newEEG = myEEG;
        newEEG.data = [];
        newEEG.data = mean(myEEG.data(j,:,:),3);
         
      %  Filter TEPs into frequency bands 
        EEG_Alpha_meanTrials = tesa_filtbutter(newEEG, 8, 12, 4, 'bandpass' );
        EEG_Beta_meanTrials = tesa_filtbutter( newEEG, 13, 30, 4, 'bandpass' );
        EEG_Gamma_meanTrials = tesa_filtbutter(newEEG, 31, 45, 4, 'bandpass' );
        EEG_BetaGamma_meanTrials = tesa_filtbutter(newEEG, 13, 45, 4, 'bandpass' );
        
        % Hilbert transformation of data at each freq band
        hilb_Alpha_meanTrials(j,:,idx) = hilbert(EEG_Alpha_meanTrials.data);
        hilb_Beta_meanTrials(j,:,idx) = hilbert(EEG_Beta_meanTrials.data);
        hilb_Gamma_meanTrials(j,:,idx) = hilbert(EEG_Gamma_meanTrials.data);
        hilb_BetaGamma_meanTrials(j,:,idx) = hilbert(EEG_BetaGamma_meanTrials.data);
        
        % Hilbert evoked amplitude
        EvokedAmp.Alpha(j,:,idx) = abs(hilb_Alpha_meanTrials(j,:,idx));
        EvokedAmp.Beta(j,:,idx) = abs(hilb_Beta_meanTrials(j,:,idx));
        EvokedAmp.Gamma(j,:,idx) = abs(hilb_Gamma_meanTrials(j,:,idx));
        EvokedAmp.BetaGamma(j,:,idx) = abs(hilb_BetaGamma_meanTrials(j,:,idx));
        
        % Evoked Amplitude_baseline removed
        EvokedAmp.Alpha_bsLineCorrected(j,:,idx) = EvokedAmp.Alpha(j,:,idx) - mean(EvokedAmp.Alpha(j,baseLine,idx),2);
        EvokedAmp.Beta_bsLineCorrected(j,:,idx) = EvokedAmp.Beta(j,:,idx) - mean(EvokedAmp.Beta(j,baseLine,idx),2);
        EvokedAmp.Gamma_bsLineCorrected(j,:,idx) = EvokedAmp.Gamma(j,:,idx) - mean(EvokedAmp.Gamma(j,baseLine,idx),2);
        EvokedAmp.BetaGamma_bsLineCorrected(j,:,idx) = EvokedAmp.BetaGamma(j,:,idx) - mean(EvokedAmp.BetaGamma(j,baseLine,idx),2);
    end
end
end

