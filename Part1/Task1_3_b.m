%% Task 1.3.b Correlation Estimation
% This task is to generate the PSD estimate of several realisations of a
% random process.
% Author: Junyu Yan
%--------------------------------------------------------------------------

clc;clear;close all;
%% Initialization
% The signal length
nSample = 2048;
% The normalized signal frequency
fs = 1;
% The sampling time
time = (0: nSample - 1) / fs;
% The number of realisation
nRe = 100;
% The frequency of sinusoidal signal
f_sine = [0.2, 0.25];
% The amplitude of sinusoidal signal
A_sin = [0.3, 0.6];
% The psd matrix of all realisation
psd_matrix = [];

%% Calculate the psd estimate and plot
subplot(1,2,1);
for iRe = 1:nRe
    % The noisy sinusoidal signal
    noisySin = A_sin(1)*sin(2*pi*f_sine(1)*time)+A_sin(2)*sin(2*pi*f_sine(2)*time)+wgn(1,nSample,0);
    [ACF_sin, lag_sin] = xcorr(noisySin,'biased');
    %Calculate the PSD
    psd_sin = real(fftshift(fft(ifftshift(ACF_sin))));
    % The normalized frequency with respect to ACF
    f_sin = lag_sin ./ (2 * nSample) * fs;
    F1=plot(f_sin,psd_sin,'c','linewidth',2);
    hold on;
    psd_matrix = [psd_matrix;psd_sin];
end
psd_mean = mean(psd_matrix);
F2=plot(f_sin,psd_mean,'b','linewidth',2);
grid on; grid minor;
set(gca,'fontsize',10);
xticks(-0.5:0.1:0.5);
title('PSD estimate (different realisation and mean)');
xlabel('Normalized frequency (\pi radians/sample)');
ylabel('Magnitude');
legend([F1,F2],'Realization','Mean');

%% The standard deviation
subplot(1,2,2);
psd_std = std(psd_matrix);
plot(f_sin,psd_std,'r','linewidth',2);
grid on; grid minor;
title('Standard deviation of the PSD estimate');
xlabel('Normalized frequency (\pi radians/sample)');
ylabel('Magnitude');
set(gca,'fontsize',10);
xticks(-0.5:0.1:0.5);