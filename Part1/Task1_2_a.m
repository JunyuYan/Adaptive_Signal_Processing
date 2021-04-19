%% Task 1.2.a Periodogram-based Methods Applied to Real-World Data
% This task is to apply one periodogram-based spectral estimation technique 
% to the sunspot time series possibly with some preprocessing of data.
% Author: Junyu Yan
%--------------------------------------------------------------------------

clc;clear;close all;
load sunspot.dat;
%% Initialize parameters
% The normalized signal frequency
fs = 1;
% The sampling time
time = sunspot(:,1);
% The original signal
Signal_org = sunspot(:,2);
% The number of samples
nSample = length(Signal_org);
% The signal after removing mean and trend
Signal_mtr = detrend(Signal_org - mean(Signal_org));
% Add a small DC value to avoid 0 in log 
Signal_addDC = Signal_org + eps;
% The signal after logarithm
Signal_log = log(Signal_addDC)-mean(log(Signal_addDC));
%% Apply periodgram-based technique 
[psd_org,f] = periodogram(Signal_org,hamming(nSample),nSample,fs);
[psd_mtr,~] = periodogram(Signal_mtr,hamming(nSample),nSample,fs);
[psd_log,~] = periodogram(Signal_log,hamming(nSample),nSample,fs);
%% Plot the signal
figure;
subplot(2,1,1);
plot(time, Signal_org,'linewidth',2);
hold on;
plot(time, Signal_mtr,'linewidth',2);
hold on;
plot(time, Signal_log,'linewidth',2);
grid on;grid minor;
title('The sunspot time series');
xlabel('Time (year)');
ylabel('Number of sunspots');
set(gca,'fontsize',10);
legend('Orignal','Mean-detrend','Log-mean')
%% Plot the power spectrum density
subplot(2,1,2);
plot(f, pow2db(psd_org),'linewidth',2);
hold on;
plot(f, pow2db(psd_mtr),'--','linewidth',2);
hold on;
plot(f, pow2db(psd_log),'linewidth',2);
grid on;grid minor;
title('The power spectrum density of sunspot series with hamming window');
xlabel('Normalised frequency (\pi radians/sample)');
ylabel('Magnitude (dB)');
set(gca,'fontsize',10);
legend('Orignal','Mean & detrend','Log & mean')
