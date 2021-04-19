%% Task 1.6.c Robust Regression
% This task is to calculate the OLS and PCR solutions, and then compare the
% error between the true value and estimated value with both training and
% testing data
% Author: Junyu Yan
%--------------------------------------------------------------------------

clc;clear;close all
load('../Data/PCAPCR.mat');

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
% Generate the low rank Xtest
[Utest_x,Stest_x,Vtest_x] = svd(Xtest);
Xtest_lowrank = Utest_x(:,1:xtest_rank)*Stest_x(1:xtest_rank,1:xtest_rank)*Vtest_x(:,1:xtest_rank)';
%% Calculating solution
% OLS method
% Train set
Ytrain_OLS =  Xnoise*B_OLS;
% Test set
Ytest_OLS = Xtest*B_OLS;
% PCR method
% Train set
Ytrain_PCR = Xnoise*B_PCR;
% Test set
Ytest_PCR = Xtest_lowrank*B_PCR;

%% Calculate the error
errorTrain_OLS = abs(vecnorm(Y-Ytrain_OLS)).^2;
errorTest_OLS = abs(vecnorm(Ytest-Ytest_OLS)).^2;
errorTrain_PCR = abs(vecnorm(Y-Ytrain_PCR)).^2;
errorTest_PCR = abs(vecnorm(Ytest-Ytest_PCR)).^2;

%% Total error
ToterrorTrain_OLS = sum(errorTrain_OLS);
display(['The total error of training data with OLS is ',num2str(ToterrorTrain_OLS)]);
ToterrorTest_OLS = sum(errorTest_OLS);
display(['The total error of testing data with OLS is ',num2str(ToterrorTest_OLS)]);
ToterrorTrain_PCR = sum(errorTrain_PCR);
display(['The total error of training data with PCR is ',num2str(ToterrorTrain_PCR)]);
ToterrorTest_PCR = sum(errorTest_PCR);
display(['The total error of testing data with PCR is ',num2str(ToterrorTest_PCR)]);

%% Plot the results
figure;
% Plot the error
subplot(1,2,1);
stem(errorTrain_OLS,'b-*');
hold on;
stem(errorTrain_PCR,'r--o');
grid on;grid minor;
xlabel('Subspace dimension index');
ylabel('Square error');
set(gca,'fontsize',10);
title('Error between reproduced and original training output');
legend('OLS','PCR');
subplot(1,2,2);
stem(errorTest_OLS,'b-*');
hold on;
stem(errorTest_PCR,'r--o');
grid on; grid minor;
xlabel('Subspace dimension index');
ylabel('Square error');
set(gca,'fontsize',10);
title('Error between reproduced and orignal testing data');
legend('OLS','PCR');
