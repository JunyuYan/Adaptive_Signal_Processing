%% Task 1.6.a Robust Regression
% This task is to obtain the singular values of X and Xnoise and identify
% the rank of the input data, then plot the square error between each
% singular value of X and Xnoise
% Author: Junyu Yan
%--------------------------------------------------------------------------

clc;clear;close all
load('PCAPCR.mat');

%% Obtain the singular values and rank
% Obtain the singular values of X
x_svd = svd(X);
% Obtain the singular values of Xnoise
xNoise_svd = svd(Xnoise);
% Obtain the rank of input signal
x_rank = rank(X);
% Obtain the square error between singular values of X and Xnoise
error = abs(x_svd-xNoise_svd).^2;

%% Plot the results
figure;
% Plot signular values
subplot(2,1,1);
stem(xNoise_svd,'b-*');
hold on;
stem(x_svd,'r--o');
grid on;grid minor;
xlabel('Subspace dimension index');
ylabel('Singular values');
title('The singular values of original input and noisy inputs');
legend('Noisy input','Original input')
% Plot the error
subplot(2,1,2);
stem(error);
grid on;grid minor;
xlabel('Subspace dimension index');
ylabel('Square error');
title('The square error between original input and noisy inputs');
legend('Square error');