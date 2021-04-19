%% Task 1.5.a Respiratory Sinus Arrhythmia from RR-Intervals
% This task is to compare the power spectrum density of the signal
% estimated by AR model with different orders and 10000 samples. And
% compare the effect of under-modeling or over-modeling
% Author: Junyu Yan
%--------------------------------------------------------------------------

clc;clear;close all
% load('RRI-DATA.mat');
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

%% Standard periodogram and averaged periodogram
for iTrial = 1:length(RRI)
    subplot(length(RRI),1,iTrial);
    %% Standard periodogarm
    % The number of samples
    nSample = length(RRI{iTrial});
    [psd_stand, f_stand] = periodogram(RRI{iTrial}, hamming(nSample), nfft, fs, 'onesided');
    % Plot the standard psd
    plot(f_stand, pow2db(psd_stand),'linewidth', 2);
    hold on;
    %% Averaged periodogram
    for index = 1:length(winLength)
        N = winLength(index) * fs;
        [psd_avg, f_avg] = pwelch(RRI{iTrial}, hamming(N), noverlap, nfft, fs);
        plot(f_avg, pow2db(psd_avg),'linewidth', 2);
        hold on;
    end
    %% label on figure
    grid on;grid minor;
    ylim([-80,0]);
    title(['The periodogram of the RRI data: Trial ',num2str(iTrial)]);
    xlabel('Frequency (Hz)');
    ylabel('Magnitude (dB)');
    set(gca,'fontsize',10);
    legend('Standard','Averaged WL=150s','Averaged WL=50s');
end  
   

