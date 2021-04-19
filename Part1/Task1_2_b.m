%% Task 1.2.b Periodogram-based Methods Applied to Real-World Data
% This task is to apply standard periodogram approach and the averaged
% periodogram with different window length to the EEG signal
% Author: Junyu Yan
%--------------------------------------------------------------------------

clc;clear;close all
load('../Data/EEG_Data_Assignment1.mat');
%% Initialization
% The number of samples
nSample = length(POz);
% The number of DFT samples, assume 5 DFT samples per Hz
nfft = fs*5;
% The window length
winLength = [10, 5, 1];
% The overlap of window
noverlap = 0;

% Remove the mean
POz = POz-mean(POz);

%% Standard periodogram
[psd_egg_stand, f_stand] = periodogram(POz, hamming(nSample), nfft, fs, 'onesided');
% Plot the standard psd
figure;
subplot(2,1,1);
plot(f_stand, pow2db(psd_egg_stand),'linewidth', 2);
grid on;grid minor;
xlim([0,60]);
ylim([-150,-90]);
title('The PSD of EEG series: standard periodogram');
xlabel('Frequency (Hz)');
ylabel('Magnitude (dB)');
set(gca,'fontsize',12);
%% Averaged periodogram
subplot(2,1,2);
for index = 1:length(winLength)
    N = winLength(index) * fs;
    [psd_egg_avg, f_avg] = pwelch(POz, hamming(N), noverlap, nfft, fs);
    plot(f_avg, pow2db(psd_egg_avg),'linewidth', 2);
    hold on;
end
grid on;grid minor;
xlim([0,60]);
ylim([-150,-90]);
title('The PSD of EEG series: averaged periodogram');
xlabel('Frequency (Hz)');
ylabel('Magnitude (dB)');
legend('Window=10s','Window=5s','Window=1s');
set(gca,'fontsize',12);

