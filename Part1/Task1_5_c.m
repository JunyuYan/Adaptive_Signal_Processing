%% Task 1.5.c Respiratory Sinus Arrhythmia from RR-Intervals
% This task is to compare the power spectrum density of the signal
% estimated by AR model with different orders and 10000 samples. And
% compare the effect of under-modeling or over-modeling
% Author: Junyu Yan
%--------------------------------------------------------------------------

clc;clear;close all
load ('../Data/xRRI1.mat');
load ('../Data/xRRI2.mat');
load ('../Data/xRRI3.mat');

%% Initialization
% RRI data
RRI = {xRRI_trial1,xRRI_trial2,xRRI_trial3};
% Remove the mean and detrend
RRI = {detrend(RRI{1}-mean(RRI{1})),detrend(RRI{2}-mean(RRI{2})),detrend(RRI{3}-mean(RRI{3}))};
% The sampling frequency
fs = 4;
% The number of DFT samples
nfft = 2048;
% The window length
winLength = [150, 50];
% The overlap of window
noverlap = 0;
% The order
order = 2:4:10;

%% Estimate the power spectrum density
for iTrial = 1:length(RRI)
    subplot(length(RRI),1,iTrial);
    %% Standard periodogarm
    % The number of samples
    nSample = length(RRI{iTrial});
    [psd_stand, f_stand] = periodogram(RRI{iTrial}, hamming(nSample), nfft, fs, 'onesided');
    % Plot the standard psd
    plot(f_stand, pow2db(psd_stand),'linewidth', 2);
    hold on;
    %% AR periodogram
    for iOrder = order
        nData = length(RRI{iTrial}); 
        % Estimate the coefficients and variance of the AR model
        [coeff_est, var_est] = aryule(RRI{iTrial}, iOrder);
        % Obtain the estimated AR model
        [H_est,w_est] = freqz(sqrt(var_est),coeff_est,nData,fs);
        plot(w_est, pow2db(abs(H_est).^2),'linewidth',2);
        hold on;
    end
    %% label on figure
    grid on;grid minor;
    ylim([-80,5]);
    title(['The periodogram of the RRI data: Trial ',num2str(iTrial)]);
    xlabel('Frequency (Hz)');
    ylabel('Magnitude (dB)');
    set(gca,'fontsize',10);
    legend('Standard','AR Order = 2','AR Order = 6','AR Order=10');
end  
   