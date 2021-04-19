%% Task 2.3.d EEG Data noise cancellation
% Author: Junyu Yan
%--------------------------------------------------------------------------

clc;clear;close all;
EEG = load('../Data/EEG_Data_Assignment2.mat');
%% Initialization
% EEG signal
eegSig = EEG.POz;
% EEG signal preprocessing
eegSig = (eegSig-mean(eegSig))';
% The number of samples of all experiments
nSample = length(eegSig);
% Sampling frequency
fs = EEG.fs;
% Learning step size
step =[0.1,0.005,0.001];
% Delay
delay = 0;
% Leakage
leakage = 0;
% Steady state offset
t0 = 50;
% Filter length
orderFilter = 5:5:20;

% The sinusoid parameters
% Time
t = (0:nSample-1)/fs;
% Amplitude
aSine = 1;
% Normalized frequency
fSine = 50;
% Sinusiod wave
xn = aSine*sin(2*pi*fSine*t);
% Noise power
varNoise = 0.01;
% Synthetic noisy input
noisySig = xn + varNoise*randn(1,nSample);

% Spectrum parameter
% The length of window
winLength = 2^12;
% The number of FFT
nFFT = 2^13;
% The percentage of overlap 
overlap = 0.8;

% Paramters initialization
% ANC predicted noise
noise_ANC = cell(length(orderFilter),length(step));
% ANC predicted signal
pred_ANC = cell(length(orderFilter),length(step));

%% ANC
for iOrder = 1:length(orderFilter)
    inputSig_ANC = noisySig;
    desireSig_ANC = [0,eegSig(1:end-1)];
    for iStep = 1:length(step)
        % Apply ANC algorithm
        [~,~,noise_ANC{iOrder,iStep}] = funLMS(inputSig_ANC,desireSig_ANC,orderFilter(iOrder),step(iStep),delay,leakage);
        % Predicted xn
        pred_ANC{iOrder,iStep} = desireSig_ANC - noise_ANC{iOrder,iStep};
    end
end

% The spectrum of original POz signal
figure;
spectrogram(eegSig,winLength,round(overlap*winLength),nFFT, fs,'yaxis'); 
ylim([0,60]);
title('Original POz spectrum');

% The spectum of POz signal after ANC
figure;
FigIndex = 1;
for iOrder = 1:length(orderFilter)-1
    for iStep = 1:length(step)
        subplot(length(step),3,FigIndex);
        FigIndex = FigIndex + 1;
        spectrogram(pred_ANC{iOrder,iStep},winLength,round(overlap*winLength),nFFT, fs,'yaxis'); 
        ylim([0,60]);
        yticks(0:10:60)
        title(sprintf('POz ANC spectrum\nM = %d, \\mu = %.3f', orderFilter(iOrder),step(iStep)));
    end
end

%% Periodogram of original and ANC POz
% The optimal filter order and step size
optimalM = 10;
optimalU = 0.001;
% Find index for optimal values
M_index = find(optimalM == orderFilter);
U_index = find(optimalU == step);
% Obtain periodgram
[psd_org, f_org] = periodogram(eegSig, hamming(nSample), nFFT, fs, 'onesided');
[psd_ANC, f_ANC] = periodogram(pred_ANC{M_index,U_index}, hamming(nSample), nFFT, fs, 'onesided');
% Plot the periodgram
figure;
subplot(1,2,1);
plot(f_org, pow2db(psd_org),'linewidth', 2);
hold on;
plot(f_ANC, pow2db(psd_ANC),'linewidth', 2);
grid on;grid minor;
xlim([0,60]);
xticks(0:10:60);
ylim([-160,-100]);
legend('Original','ANC');
title('Periodgram of original and ANC POz');
xlabel('Frequency (Hz)');
ylabel('Magnitude (dB)');
% Plot the periodgram error
subplot(1,2,2);
plot(f_org, pow2db(abs(psd_org-psd_ANC)),'linewidth', 2);
grid on;grid minor;
xlim([0,60]);
xticks(0:10:60);
ylim([-160,-100]);
title('Periodgram error between original and ANC POz');
xlabel('Frequency (Hz)');
ylabel('Magnitude (dB)');