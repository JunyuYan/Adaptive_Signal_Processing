%% Task 1.6.c Robust Regression
% This task is to compare the mean square error estimates of OLS and PCR 
% with regval.
% Author: Junyu Yan
%--------------------------------------------------------------------------

clc;clear;close all
load('../Data/PCAPCR.mat');

%% Initialization
% The number of iteration
N_iter = 100;
% The OLS error
error_OLS = zeros(N_iter,size(Y,2));
% The PCR error
error_PCR = zeros(N_iter,size(Y,2));

%% Calculate regression matrix B 
% OLS method
B_OLS = (Xnoise'*Xnoise)\Xnoise'*Y;
% PCR method
% Calculate the rank
x_rank = rank(X);
xtest_rank = rank(Xtest);
% Generate the matrix B 
[U_xNoise,S_xNoise,V_xNoise] = svd(Xnoise);
B_PCR = V_xNoise(:,1:x_rank)/S_xNoise(1:x_rank,1:x_rank)*U_xNoise(:,1:x_rank)'*Y;

%% Calculate the error
for i_iter = 1 : N_iter
    [Yest_OLS,Ytest_OLS] = regval(B_OLS);
    error_OLS(i_iter,:) = abs(vecnorm(Ytest_OLS-Yest_OLS)).^2;
    [Yest_PCR,Ytest_PCR] = regval(B_PCR);
    error_PCR(i_iter,:) = abs(vecnorm(Ytest_PCR-Yest_PCR)).^2;
end
MSE_OLS = mean(error_OLS);
MSE_PCR = mean(error_PCR);
% Total error
Toterror_OLS = sum(MSE_OLS);
display(['The total error with OLS is ',num2str(Toterror_OLS)]);
Toterror_PCR = sum(MSE_PCR);
display(['The total error with PCR is ',num2str(Toterror_PCR)]);

%% Plot the results
figure;
stem(MSE_OLS,'b-*');
hold on;
stem(MSE_PCR,'r--o');
grid on;grid minor;
xlabel('Subspace dimension index');
ylabel('Square error');
title('Error between reproduced test data and its estimation with OLS and PCR coefficients');
legend('OLS','PCR');

