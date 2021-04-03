%% Task 3.2.d Time frequency spectrum for EEG signal by DFT-CLMS
% Author: Junyu Yan
%--------------------------------------------------------------------------

clc;clear;close all;
EEG = load('EEG_Data_Assignment1.mat');
%% Initialization
% EEG signal
EEG_Sig = EEG.POz;
% Segment start
Sig_init = 1000;
% The number of samples
nSample = 1200;
% Proessed EEG 
EEG_Sig = EEG_Sig(Sig_init:Sig_init+nSample-1)';
% Remove mean and detrend signal
EEG_Sig = detrend(EEG_Sig-mean(EEG_Sig));
% The sampling frequency
fs = EEG.fs;
% Order
orderAR = 1;
% Frequency for time-frequency diagram
freq = (0:nSample-1).*(fs/nSample);
% CLMS parameters
% Learning step size
step = 1;
% Leakage
leakage = [0,0.001];
% Delays
delay = 1;
%% CLMS
figure;
for iLeak= 1:length(leakage)
    inputSig = (1/nSample).*exp(1i*2*pi*(0:nSample-1).'*(1:nSample)/nSample);
    desireSig = EEG_Sig;
    % Apply CLMS algorithm
    [weight_CLMS,~,~] = funDFT_CLMS(inputSig,desireSig,step,leakage(iLeak));
    H = abs(weight_CLMS).^2;
    % Remove outliers
    medianH = 1000 * median(median(H));
    H(H>medianH) = medianH;
    %% Plot the time-frequency figure
    subplot(2,1,iLeak);
    surf(1:nSample,freq,H,'LineStyle','none');
    view(2);
    cbar = colorbar;
    cbar.Label.String = 'PSD(dB)';
    title(sprintf('EEG signal spectrogram, AR(%d),\\gamma=%.3f',orderAR,leakage(iLeak)));
    xlabel('Time (Sample)');
    ylabel('Frequency (Hz)');
    ylim([0,70]);
end