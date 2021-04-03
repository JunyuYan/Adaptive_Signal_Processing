%% Task 3.2.a The power spectrum for frequency modulated (FM) signal
% Author: Junyu Yan
%--------------------------------------------------------------------------

clc;clear;close all;
%% Initialization
% The number of samples
nSample = 1500;
% The sampling frequency
fs = 1500;
% The time 
t = (0:nSample-1)/fs;
% Noise power
varNoise = 0.05;
% Order
orderAR = [1,5,10];
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

%% Plot frequency and phase of FM signal
% The frequency of FM signal
figure;
subplot(1,2,1);
plot(fn,'LineWidth',2);
grid on;grid minor;
title('Frequency of FM signal');
xlabel('Time (Sample)');
ylabel('Frequency (dB)');
set(gca,'fontsize',12);
% The phase of FM signal
subplot(1,2,2);
plot(angle(exp(1i*(2*pi*phase/fs))),'LineWidth',2);
grid on;grid minor;
title('Phase of FM signal');
xlabel('Frequency (Hz)');
ylabel('Phase (rad)');
set(gca,'fontsize',12);


%% Power spectrum for FM signal
figure;
for iOrder = 1:length(orderAR)
    % Estimate AR coefficient
    a_coeff = aryule(FM_Sig,orderAR(iOrder));
    % Frequency response
    [h_freq,w_freq] = freqz(1,a_coeff,nSample,fs);
    % Power spectrum for FM signal
    psd_FM = abs(h_freq.^2);
    subplot(1,length(orderAR),iOrder);
    plot(w_freq,pow2db(psd_FM),'LineWidth',2);
    grid on; grid minor;
    title(sprintf('Power spectrum for FM signal with AR(%d)',orderAR(iOrder)));
    xlabel('Frequency (Hz)');
    ylabel('PSD (dB)');
    set(gca,'fontsize',12);
    xticks(0:100:800);
end
%% Segemented power spectrum for FM signal
% Segment
Seg = 500;
% Number of segments
nSeg = 3;
figure;
for iSeg = 1:nSeg
    Seg_Sig = FM_Sig((iSeg-1)*Seg+1:iSeg*Seg);
    % Coefficients for segmented signal
    a_coeff = aryule(Seg_Sig,1);
    % Frequency response
    [h_freq,w_freq] = freqz(1,a_coeff,nSample/nSeg,fs);
    % Power spectrum for FM signal
    psd_FM = abs(h_freq.^2);
    subplot(1,nSeg,iSeg);
    plot(w_freq,pow2db(psd_FM),'LineWidth',2);
    grid on; grid minor;
    title(sprintf('FM segment %d: Power spectrum AR(1)',iSeg));
    xlabel('Time(Sample)');
    ylabel('PSD (dB)');
    set(gca,'fontsize',12);
    xticks(0:100:800);
end
