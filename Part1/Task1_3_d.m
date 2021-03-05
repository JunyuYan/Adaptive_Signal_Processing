%% Task 1.3.d Correlation Estimation
% This task is to generate signals composed of exponential signals with
% different length and frequencies, then obtain the psd to verify that the
% periodogram will show correct line spectral by considering more data
% samples
% Author: Junyu Yan
%--------------------------------------------------------------------------

clc;clear;close all;
%% Initialization
% The number of samples
N = [20,30,40,50];
% The normalized signal frequency
fs = 1;
% The frequency of sinusoidal signal
f_sin = [0.3, 0.32, 0.35];
% The amplitude of sinusoidal signal
A_sin = [1, 1, 1];
% Noise power
nPower = 0.2;
% The number of dft
nDFT = 512;

figure;
%% Periodogarm for signal with one exponential
legendstr = [];
subplot(3,1,1);
for i = 1:length(N)
    % The sampling time
    n = (0: N(i) - 1) / fs;
    % Generate noise
    noise = sqrt(nPower/2)*(randn(size(n))+1j*randn(size(n)));
    % Generate signal
    x1 = A_sin(1)*exp(1j*2*pi*f_sin(1)*n)+noise;
    [psd1,f1] = periodogram(x1,rectwin(N(i)),nDFT, fs);
    plot(f1,psd1,'linewidth',2);
    legendstr = [legendstr; 'N=',num2str(N(i))];
    hold on;
end
grid on;grid minor;
xlim([0.25,0.4]);
xlabel('Frequency (Hz)');
ylabel('PSD (dB)');
title('Periodogram Power Spectral Density Estimate for one exponentials')
legend(legendstr)

%% Periodogram for signal with two exponentials
legendstr = [];
subplot(3,1,2);
for i = 1:length(N)
    % The sampling time
    n = (0: N(i) - 1) / fs;
    % Generate noise
    noise = sqrt(nPower/2)*(randn(size(n))+1j*randn(size(n)));
    %Generate signal
    x2 = A_sin(1)*exp(1j*2*pi*f_sin(1)*n)+A_sin(2)*exp(1j*2*pi*f_sin(2)*n)+noise;
    [psd2,f2] = periodogram(x2,rectwin(N(i)),nDFT, fs);
    plot(f2,psd2,'linewidth',2);
    legendstr = [legendstr; 'N=',num2str(N(i))];
    hold on;
end
grid on; grid minor;
xlim([0.25,0.4]);
xlabel('Frequency (Hz)');
ylabel('PSD (dB)');
title('Periodogram Power Spectral Density Estimate for two exponentials');
legend(legendstr);

%% Periodogram for signal with three exponentials
legendstr = [];
subplot(3,1,3);
for i = 1:length(N)
    % The sampling time
    n = (0: N(i) - 1) / fs;
    % Generate noise
    noise = sqrt(nPower/2)*(randn(size(n))+1j*randn(size(n)));
    % Generate signal
    x3 = A_sin(1)*exp(1j*2*pi*f_sin(1)*n)+A_sin(2)*exp(1j*2*pi*f_sin(2)*n)+A_sin(3)*exp(1j*2*pi*f_sin(3)*n)+noise;
    [psd3,f3] = periodogram(x3,rectwin(N(i)),nDFT, fs);
    plot(f3,psd3,'linewidth',2);
    legendstr = [legendstr; 'N=',num2str(N(i))];
    hold on;
end
grid on;grid minor;
xlim([0.25,0.4]);
xlabel('Frequency (Hz)');
ylabel('PSD (dB)');
title('Periodogram Power Spectral Density Estimate for three exponentials');
legend(legendstr);