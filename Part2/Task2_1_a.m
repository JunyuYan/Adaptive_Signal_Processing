%% Task 2.1.b LMS Algorithm
% This task is to implement an LMS adaptive predictor and plot the squared
% prediction error and the learning curve
% Author: Junyu Yan
%--------------------------------------------------------------------------

clc;clear;close all
%% Initialization
% The number of samples of all experiments
nSample = 1000;
% The number of realizations
nReal = 100;
% The step size of all experiments
step = [0.05 ; 0.01];
% Delay
delay = 1;
% Leakage coefficient
leakage = 0;
% The AR process parameters
% AR process coefficents
AR_coeff = [0.1, 0.8];
% The AR process order
order = length(AR_coeff);
% The noise power
varNoise = 0.25;
% The error
error = zeros(nReal,nSample,length(step));
% Average error 
av_error = zeros(length(step),nSample);
%% Generate the AR process
AR_model = arima('AR',AR_coeff,'Variance',varNoise,'Constant',0);
xn = simulate(AR_model,nSample,'NumPaths',nReal).';

%% LMS
for iStep = 1:length(step)
    for iReal = 1:nReal
        inputSig = xn(iReal,:);
        desireSig = xn(iReal,:);
        [~,error(iReal,:,iStep),~] = funLMS(inputSig,desireSig,order,step(iStep),delay,leakage);
    end
    av_error(iStep,:) = mean(error(:,:,iStep).^2);
end

%% Plot the result
figure;
subplot(1,2,1);
for iStep = 1:length(step)
    plot(pow2db(error(1,:,iStep).^2),'LineWidth',2);
    hold on;
end
xlabel('Time (Sample)');
ylabel('The squared prediction error (dB)');
title('The squared error by adaptive LMS: one realisation');
set(gca,'fontsize',12);
grid on; grid minor;
legend('step = 0.05','step = 0.01');
set(gca,'fontsize',12);
subplot(1,2,2);
for iStep = 1:length(step)
    plot(pow2db(av_error(iStep,:)),'LineWidth',2);
    hold on;
end
xlabel('Time (Sample)');
ylabel('The squared prediction error (dB)');
title('The learning curve of adaptive LMS');
grid on; grid minor;
legend('step = 0.05','step = 0.01');
set(gca,'FontSize',12);