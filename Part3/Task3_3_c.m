%% Task 3.2.c Time frequency spectrum for (FM) signal by DFT-CLMS
% Author: Junyu Yan
%--------------------------------------------------------------------------

clc;clear;close all;
%% Initialization
% The number of samples
nSample = 1500;
% The sampling frequency
fs = 1500;
% Noise power
varNoise = 0.05;
% Order
orderAR = 1;
% The circular complex valued white noise
noise = sqrt(varNoise/2)*(randn(1,nSample)+1i*randn(1,nSample));
% Function to determine the frequency of FM signal
funFreq = @(n)(((1<=n)&(n<=500)).*100+((501<=n)&(n<=1000)).*(100+(n-500)/2)+((1001<=n)&(n<=1500)).*(100+((n-1000)/25).^2));
% The frequency
fn = funFreq(1:nSample);
% The phase
phase = cumsum(fn);
% The FM signal
FM_Sig = exp(1i*(2*pi*phase/fs))+noise;
% Frequency for time-frequency diagram
freq = (0:nSample-1).*(fs/nSample);

% CLMS parameters
% Learning step size
step = 1;
% Leakage
leakage = [0,0.01,0.1,0.5];
% Delays
delay = 1;
%% CLMS
figure;
for iLeak= 1:length(leakage)
    inputSig = (1/nSample).*exp(1i*2*pi*(0:nSample-1).'*(1:nSample)/nSample);
        desireSig = FM_Sig;
    % Apply CLMS algorithm
    [weight_CLMS,~,~] = funDFT_CLMS(inputSig,desireSig,step,leakage(iLeak));
    H = abs(weight_CLMS).^2;
    % Remove outliers
    medianH = 50 * median(median(H));
    H(H>medianH) = medianH;
    %% Plot the time-frequency figure
    subplot(2,2,iLeak);
    surf(1:nSample,freq,H,'LineStyle','none');
    view(2);
    cbar = colorbar;
    cbar.Label.String = 'PSD(dB)';
    title(sprintf('FM signal spectrogram, AR(%d),\\gamma=%.2f',orderAR,leakage(iLeak)));
    xlabel('Time (Sample)');
    ylabel('Frequency (Hz)');
    ylim([0,1000]);
end

%% Zero padding K-point DFT-CLMS
% Number of k points
K = 2048;
% Frequency for time-frequency diagram
freq = (0:K-1).*(fs/K);
figure;
inputSig = (1/nSample).*exp(1i*2*pi*(0:K-1).'*(1:nSample)/nSample);
desireSig = FM_Sig;
% Apply CLMS algorithm
[weight_CLMS,~,~] = funDFT_CLMS(inputSig,desireSig,step,leakage(1));
H = abs(weight_CLMS).^2;
% Remove outliers
medianH = 50 * median(median(H));
H(H>medianH) = medianH;
%% Plot the time-frequency figure
surf(1:nSample,freq,H,'LineStyle','none');
view(2);
cbar = colorbar;
cbar.Label.String = 'PSD(dB)';
title(sprintf('FM signal spectrogram, AR(%d),\\gamma=%.2f',orderAR,leakage(1)));
xlabel('Time (Sample)');
ylabel('Frequency (Hz)');
ylim([0,1000]);

