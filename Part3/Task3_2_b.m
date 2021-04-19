%% Task 3.2.b Time frequency spectrum for (FM) signal by CLMS
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
% CLMS parameters
% Learning step size
step = [0.001,0.01,0.1,1];
% Leakage
leakage = 0;
% Delays
delay = 1;

%% CLMS
figure;
for iStep= 1:length(step)
    inputSig = FM_Sig;
    desireSig = FM_Sig;
    % Apply CLMS algorithm
    [weight_CLMS,~,~] = funCLMS(inputSig,desireSig,orderAR,step(iStep),delay,leakage);
    % Power spectrum
    H = zeros(1024,nSample);
    for n = 1:nSample
        % Compute power spectrum
        [h,w] = freqz(1,[1;-conj(weight_CLMS(n))],1024,fs);
        % Store it in a matrix
        H(:,n) = abs(h).^2;
    end
    % Remove outliers
    medianH = 50 * median(median(H));
    H(H>medianH) = medianH;
    %% Plot the time-frequency figure
    subplot(2,2,iStep);
    surf(1:nSample,w,H,'LineStyle','none');
    view(2);
    cbar = colorbar;
    cbar.Label.String = 'PSD(dB)';
    title(sprintf('FM signal spectrogram, AR(%d),\\mu=%.3f',orderAR,step(iStep)));
    xlabel('Time (Sample)');
    ylabel('Frequency (Hz)');
    ylim([0,600]);
    set(gca,'fontsize',10);
end
