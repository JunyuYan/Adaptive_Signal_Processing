%% Task 1.4.b Spectrum of Autoregressive Processes
% This task is to compare the power spectrum density of the signal
% estimated by AR model with different orders and 1000 samples
% Author: Junyu Yan
%--------------------------------------------------------------------------

clc; clear; close all;
%% Initialize parameters
% The number of samples
nSample = 1000;
% The coefficients for AR model
coeffAR = [2.76, -3.81, 2.65, -0.92];
% The variance for AR model
var = 1;
% The constant for AR model
const = 0;
% The number of sample to discard
nDiscard = 500;
% Build the data AR model
AR = arima('Constant',const,'AR',coeffAR,'Variance',var);
% Generate data
data = simulate(AR,nSample);
% Discard the first 500 samples
data = data(nDiscard:end);
% The new data length
nData = length(data);
% The ideal AR model
[H,w] = freqz(sqrt(var),[1,-coeffAR],nData);
% The order of AR model
order = 2:14;
% The noise power
nPower = zeros(1,length(order));
% The power spectrum matrix
psd_matrix = zeros(length(data),length(order));
%% Obtain AR estimated PSD
for iOrder = order
    % Estimate the coefficients and variance of the AR model
    [coeff_est, var_est] = aryule(data, iOrder);
    % Obtain the estimated AR model
    [H_est,w_est] = freqz(sqrt(var_est),coeff_est,nData);
    % Store the psd
    psd_matrix(:,iOrder) = pow2db(abs(H_est.^2));
    % Store the power of noise signal
    nPower(:,iOrder-1) = var_est; 
end
subplot(2,1,1);
plot(order,nPower,'linewidth',2);
grid on;grid minor;
xlabel('Normalised frequency (\pi radians/sample)');
ylabel('Magnitude (dB)');
title('The Prediction error with different orders');
subplot(2,1,2);
%% Plot the PSD 
% The legend string
plot(w,pow2db(abs(H).^2),'linewidth',2);
hold on;
for iOrder = 6:8
    plot(w_est,psd_matrix(:,iOrder),'linewidth',2);
    hold on;
end
grid on;grid minor;
xlabel('Normalised frequency (\pi radians/sample)');
ylabel('PSD (dB)');
title('The Power Spectrum Density with different orders and 950 samples');
legend('Ground Truth','Order=6','Order=7','Order=8');



