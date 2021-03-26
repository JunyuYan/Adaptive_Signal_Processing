%% Task 2.1.c LMS misadjustment
% This task is to estimate the LMS midadjustment by steady state mean 
% square error and compare it with the theorectical one 
% Author: Junyu Yan
%--------------------------------------------------------------------------

clc;clear;close all
%% Initialization
% Steady state time index
steady_t0 = 500;
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
% The auto-correlation matrix of input
corr = [25/27, 25/54; 25/54, 25/27];
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
% The MSE error
MSE_error = zeros(length(step),nReal);

%% Generate the AR process
AR_model = arima('AR',AR_coeff,'Variance',varNoise,'Constant',0);
xn = simulate(AR_model,nSample,'NumPaths',nReal).';

%% LMS
for iStep = 1:length(step)
    for iReal = 1:nReal
        inputSig = xn(iReal,:);
        desireSig = xn(iReal,:);
        [~,error(iReal,:,iStep),~] = funLMS(inputSig,desireSig,order,step(iStep),delay,leakage);
        MSE_error(iStep,iReal) = mean(error(iReal,steady_t0:end,iStep).^2);
    end
end

% Obtain EMSE error
EMSE_error = mean(MSE_error-varNoise,2);
% Obtain the estimated LMS misadjustment
M_est = EMSE_error/varNoise;

% Obtain the theoretical LMS misadjustment
M_theo = step./2 * trace(corr);

%% Print the result
for iStep = 1:length(step)
    % Print the theoretical adjustment
    fprintf("[step size: %.2f], theoretical misadjustment: %.4f\n", step(iStep), M_theo(iStep));
    % Print the estimated adjustment
    fprintf("[step size: %.2f], estimated misadjustment: %.4f\n", step(iStep), M_est(iStep));
end