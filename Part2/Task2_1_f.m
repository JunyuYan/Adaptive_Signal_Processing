%% Task 2.1.f LMS steady state values of the adaptive filter coefficients
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
leakage = 0.2:0.3:0.8;
% The AR process parameters
% AR process coefficents
AR_coeff = [0.1, 0.8];
% The AR process order
order = length(AR_coeff);
% The noise power
varNoise = 0.25;
% The weight
weight = cell(length(step),nReal);
% Average error 
av_weight = cell(length(step),1);

%% Generate the AR process
AR_model = arima('AR',AR_coeff,'Variance',varNoise,'Constant',0);
xn = simulate(AR_model,nSample,'NumPaths',nReal).';

%% LMS
for iLeak = 1:length(leakage)
    for iStep = 1:length(step)
        for iReal = 1:nReal
            inputSig = xn(iReal,:);
            desireSig = xn(iReal,:);
            [weight{iStep,iReal},~,~] = funLMS(inputSig,desireSig,order,step(iStep),delay,leakage(iLeak));
        end
        av_weight{iStep} = mean(cat(3,weight{iStep,:}),3);

        %% Plot the result
        subplot(length(leakage),length(step),2*(iLeak-1)+iStep);
        for iOrder = 1:order
            plot(av_weight{iStep}(iOrder,:),'LineWidth',2);
            hold on;
            plot([1, nSample],[AR_coeff(iOrder) AR_coeff(iOrder)],'--','LineWidth',2);
        end
        title(sprintf('Leakly LMS Weights: \\gamma= %.1f, \\mu = %.2f', leakage(iLeak),step(iStep)));
        ylim([-0.2 1.2]);
        xlabel('Number of iterations (Samples)');
        ylabel('Average weight');
        legend('a1\_Est','a1','a2\_Est','a2','Location','North','NumColumns',2*order);
        grid on; grid minor;
        set(gca,'fontsize',12);
    end
end