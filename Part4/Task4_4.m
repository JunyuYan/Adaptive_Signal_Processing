%% Task 4.4 LMS algorithm with bias tanh for time-series data
% This task is to implement an LMS adaptive predictor with bias
% tanh activation function for time series data and plot the squared 
% prediction error and the prediction gain
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
% The scaling factor
scale = 10:10:100;
% The bias
bias = 1;
% Prediction error
error = zeros(length(scale),nSample);
% Mean squared error
av_error = zeros(length(scale),1);
% Prediction output signal
pred = zeros(length(scale),nSample);
% Prediction gain
Rp_gain = zeros(length(scale),1);

%% LMS
inputSig = Signal;
desireSig = Signal;
for iScale = 1:length(scale)
    [~,error(iScale,:),pred(iScale,:)] = funPerception(inputSig,desireSig,order,step,delay,leakage,scale(iScale),bias);
    % Mean-square error
    av_error(iScale)= mean(abs(error(iScale,:)).^2);
    % Prediction gain
    Rp_gain(iScale) = pow2db(var(pred(iScale,:))/var(error(iScale,:)));
end
%% Plot the result
% Plot the 
figure;
iFig = 1;
for iScale = 1:3:10
    subplot(4,1,iFig);
    iFig = iFig+1;
    plot(Signal,'k');
    hold on;
    plot(pred(iScale,:),'r');
    xlabel('Time (Sample)');
    ylabel('Ampitude');
    title(sprintf('Original signal and output signal by biased & %d scaled tanh-LMS',scale(iScale)));
    legend('Zero mean','tanh-LMS');
    set(gca,'fontsize',12);
    grid on; grid minor;
end
% Plot the mean squred error and prediction gain
figure;
yyaxis left;
plot(scale,pow2db(av_error),'LineWidth',2);
ylabel('Mean Square Error (dB)');
yyaxis right;
plot(scale,Rp_gain,'LineWidth',2);
title('The MSE and R_{p} gain for scaled & biased tanh-LMS');
xlabel('Activation scale');
ylabel('Prediction gain (dB)');
set(gca,'xtick',(10:10:100));
grid on;
grid minor;
legend('MSE','Gain');
set(gca,'fontsize',12);