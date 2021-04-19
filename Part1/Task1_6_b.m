%% Task 1.6.b Robust Regression
% This task is to use the rank of noiseless inputs to generate a low rank
% noisy inputs, and compare the difference between the noiseless inputs
% with noisy inputs and low rank noisy inputs
% Author: Junyu Yan
%--------------------------------------------------------------------------

clc;clear;close all
load('../Data/PCAPCR.mat');

%% Obtain the singular values and rank
% Obtain the rank of input signal
x_rank = rank(X);
% Obtain the singular values of Xnoise
[U_xNoise,S_xNoise,V_xNoise] = svd(Xnoise);
% Reconstruct the low rank Xnoise
Xnoise_lowrank = U_xNoise(:,1:x_rank)*S_xNoise(1:x_rank,1:x_rank)*V_xNoise(:,1:x_rank)';

%% Obtain the error between variables (columns) of X and Xnoise, X and low rank Xnoise
error_xN = abs(vecnorm(X-Xnoise)).^2;
error_xNlowrank = abs(vecnorm(X-Xnoise_lowrank)).^2;

%% Plot the results
figure;
% Plot the error
stem(error_xN,'b-*');
hold on;
stem(error_xNlowrank,'r--o');
grid on;grid minor;
xlabel('Subspace dimension index');
ylabel('Square error');
set(gca,'fontsize',10);
title('Error between noiseless inputs with noisy inputs and low rank noisy inputs');
legend('Noisy inputs','Low rank noisy inputs');
