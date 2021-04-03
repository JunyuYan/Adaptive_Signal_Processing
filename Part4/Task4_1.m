%% Task 4.1 LMS algorithm for time-series data
% This task is to implement an LMS adaptive predictor for time series data
% and plot the squared prediction error and the prediction gain
% Author: Junyu Yan
%--------------------------------------------------------------------------

clc;clear;close all
load('time-series.mat');
%% Initialization
% The time series data
Signal = y';
% The number of samples of all experiments
nSample = length(Signal);
% The learning rate 
step = 1e-5;
% Delay
delay = 1;
% Leakage coefficient
leakage = 0;
% The AR process order
order = 4;

%% Preprocessing data
% Remove main
Signal = Signal-mean(Signal);

%% LMS
inputSig = Signal;
desireSig = Signal;
[~,error,pred] = funLMS(inputSig,desireSig,order,step,delay,leakage);
% Mean-square error
av_error= mean(abs(error).^2);
% Prediction gain
Rp_gain = pow2db(var(pred)/var(error));
%% Plot the result
figure;
plot(Signal,'k');
hold on;
plot(pred,'r');
xlabel('Time (Sample)');
ylabel('Ampitude');
title('Zero mean signal and output signal predicted by LMS');
legend('Zero mean','Output');
set(gca,'fontsize',12);
grid on; grid minor;
% Print the mean squred error and prediction gain
fprintf('MSE: %.4f dB \n',pow2db(av_error));
fprintf('Prediction gain: %.4f dB \n', Rp_gain);
