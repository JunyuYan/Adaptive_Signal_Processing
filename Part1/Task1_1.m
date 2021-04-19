%% Task1.1 Properties of Power Spectral Density (PSD)
% This task is to implement the approximation in the definition of PSD,
% which means compare the two definitions of the two definitions of PSD,
% and conclude when the equivalence holds or dose not hold.
% Author: Junyu Yan
% -----------------------------------------------------------------------

clc; clear; close all;
%% Initialize parameters
% The number of samples
nSample = 2048;
% The normalized s frequency
fs = 1;
% The sampling time
time = (0: nSample - 1) / fs;
% Simulation the equivalance
% For the equvalence holds, the covaranice matrix must decay
% rapidly,therefore, the impulse signal is chosen for simulation
signal_A = [1, zeros(1, nSample-1)];
% For the equvalence does not hold, the covaranice matrix decay slowly,
% therefore, the sine signal is chosen
% The normalized frequency for sine signal
f_sine = [0.5, 0.7];
signal_B = sin(2 * pi * f_sine(1) * time) + sin(2 * pi * f_sine(2) * time);
%% For definition 1: DTFT of the ACF
% Calculate the autocovaranice function
[ACF_A,lag_A] = xcorr(signal_A,'biased');
[ACF_B,lag_B] = xcorr(signal_B,'biased');
% The normalized frequency with respect to ACF
fA1 = lag_A./ (2 * nSample) * fs;
fB1 = lag_B./(2 * nSample) * fs;
% Calculate the PSD, and then shift to the center of the frequency
psd_A1 = abs(fftshift(fft(ACF_A)));
psd_B1 = abs(fftshift(fft(ACF_B)));
%% For definition 2: distribution of the average signal power over
% frequencies
% Calculate power spectral density
psd_A2 = abs(fftshift(fft(signal_A))).^2 / nSample;
psd_B2 = abs(fftshift(fft(signal_B))).^2 / nSample;
% The fequency with symmetry at center
fA2 = (-nSample/2 : nSample/2 - 1) / nSample * fs;
fB2 = (-nSample/2 : nSample/2 - 1) / nSample * fs;
%% Plot the ACF of two signals
figure;
subplot(2,1,1);
plot(lag_A, ACF_A,'LineWidth',2);
title('The autocovariance function of the impulse signal');
xlabel('Lags(sample)');
ylabel('ACF');
grid on;grid minor;
set(gca,'fontsize',10);
subplot(2,1,2)
plot(lag_A, ACF_B,'LineWidth',2);
title('The autocovariance function of the summed sine signal');
xlabel('Lags(sample)');
ylabel('ACF');
grid on;grid minor;
set(gca,'fontsize',10);
%% Plot the PSD of impulse signal for equvalance simulation
figure;
subplot(2,1,1);
plot(fA1,pow2db(psd_A1),'LineWidth',2);
hold on;
plot(fA2,pow2db(psd_A2),'--','LineWidth',2);
grid on;grid minor;
legend('Definition 1', 'Definition 2');
xlabel('Normalized frequency (\pi rad/sample)')
ylabel('Magnitude (dB)');
set(gca,'fontsize',10);
title('The power density function for two definitions of impulse signal');
%% Plot the PSD of sinusoidal signal for not equvalance simulation
subplot(2,1,2);
plot(fB1,pow2db(psd_B1),'LineWidth',2);
hold on;
plot(fB2,pow2db(psd_B2),'LineWidth',2);
grid on;
legend('Definition 1', 'Definition 2');
xlabel('Normalized frequency (\pi rad/sample)')
ylabel('Magnitude (dB)');
set(gca,'fontsize',10);
title('The power density function for two definitions of sinusoid signal');




