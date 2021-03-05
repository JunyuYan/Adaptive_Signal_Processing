%% Task 1.3.a Correlation Estimation
% This task is to calculate the baised and unbaised ACF estimate of a 
% signal, then use these ACF to compute the correlogram. Also compute the
% PSD with different ACF
% Author: Junyu Yan
%--------------------------------------------------------------------------

clc;clear;close all;
%% Initialization
% The number of samples
nSample = 2048;
% The normalized signal frequency
fs = 1;
% The sampling time
time = (0: nSample - 1) / fs;
% The frequency of sinusoidal signal
f_sine = [0.5, 0.7];
% The filter size for filtering WGN
filterSize = 5;

%% The definition of three kinds of signals
% The WGN signal with unit power
wgnSignal = wgn(1, nSample, 0);
% The noisy sinusoidal signal
noisySin = sin(2*pi*f_sine(1)*time)+sin(2*pi*f_sine(2)*time)+wgn(1,nSample,0);
% The filtered WGN signal
% Use the Moving average filter
b = 1/filterSize * ones(1,filterSize);
a = 1;
wgnFilter = filter(b,a,wgnSignal);

%% Biased and unbiased ACF
% The biased and unbiased ACF estimate of wgn signal
[ACF_wgnBi, lagBiased_wgn] = xcorr(wgnSignal,'biased');
[ACF_wgnUnBi, ~] = xcorr(wgnSignal,'unbiased');
% The biased and unbiased ACF estimate of noisey sinusoidal signal
[ACF_sinBi, lagBiased_sin] = xcorr(noisySin,'biased');
[ACF_sinUnBi, ~] = xcorr(noisySin,'unbiased');
% The biased and unbiased ACF estimate of filtered wgn signal
[ACF_filBi, lagBiased_fil] = xcorr(wgnFilter,'biased');
[ACF_filUnBi, ~] = xcorr(wgnFilter,'unbiased');

%% Plot the ACF
% Plot the baised and unbaised ACF for three signals
figure;
subplot(3,1,1);
plot(lagBiased_wgn,ACF_wgnUnBi,'linewidth',2);
hold on;
plot(lagBiased_wgn,ACF_wgnBi,'linewidth',2);
title('Correlogram of wgn signal');
xlabel('Lags (sample)');
ylabel('ACF');
grid on; grid minor;
legend('Unbiased','Biased');
subplot(3,1,2);
plot(lagBiased_sin,ACF_sinUnBi,'linewidth',2);
hold on;
plot(lagBiased_sin,ACF_sinBi,'linewidth',2);
title('Correlogram of noisy sinusoidal signal');
xlabel('Lags (sample)');
ylabel('ACF');
grid on; grid minor;
legend('Unbiased','Biased');
subplot(3,1,3);
plot(lagBiased_fil,ACF_filUnBi,'linewidth',2);
hold on;
plot(lagBiased_fil,ACF_filBi,'linewidth',2);
title('Correlogram of filtered wgn signal');
xlabel('Lags (sample)');
ylabel('ACF');
grid on; grid minor;
legend('Unbiased','Biased');

%% Obtain and plot the PSD with biased and unbiased ACF
% The normalized frequency with respect to ACF
f_wgn = lagBiased_wgn ./ (2 * nSample) * fs;
f_sin = lagBiased_sin ./ (2 * nSample) * fs;
f_fil = lagBiased_fil ./ (2 * nSample) * fs;
% Calculate the PSD
% Shift to the original frequency, then apply fft, then shift to zero 
% frequency at the center
psdBi_wgn = real(fftshift(fft(ifftshift(ACF_wgnBi))));
psdUnBi_wgn = real(fftshift(fft(ifftshift(ACF_wgnUnBi))));
psdBi_sin = real(fftshift(fft(ifftshift(ACF_sinBi))));
psdUnBi_sin = real(fftshift(fft(ifftshift(ACF_sinUnBi))));
psdBi_fil = real(fftshift(fft(ifftshift(ACF_filBi))));
psdUnBi_fil = real(fftshift(fft(ifftshift(ACF_filUnBi))));
figure;
subplot(3,1,1);
plot(f_wgn,psdUnBi_wgn,'linewidth',2);
hold on;
plot(f_wgn,psdBi_wgn,'linewidth',2);
grid on; grid minor;
title('The PSD of wgn signal with biased and unbiased ACF');
xlabel('Normalized frequency (\pi radians/sample)');
ylabel('Magnitude');
legend('Unbiased','Biased');
subplot(3,1,2);
plot(f_sin,psdUnBi_sin,'linewidth',2);
hold on;
plot(f_sin,psdBi_sin,'linewidth',2);
grid on; grid minor;
title('The PSD of noisy sinusoidal signal with biased and unbiased ACF');
xlabel('Normalized frequency (\pi radians/sample)');
ylabel('Magnitude');
legend('Unbiased','Biased');
subplot(3,1,3);
plot(f_fil,psdUnBi_fil,'linewidth',2);
hold on;
plot(f_fil,psdBi_fil,'linewidth',2);
grid on; grid minor;
title('The PSD of filtered wgn signal with biased and unbiased ACF');
xlabel('Normalized frequency (\pi radians/sample)');
ylabel('Magnitude');
legend('Unbiased','Biased');

