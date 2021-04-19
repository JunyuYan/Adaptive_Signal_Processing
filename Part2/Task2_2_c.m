%% Task 2.2.c GNGD algorithms
% Author: Junyu Yan
%--------------------------------------------------------------------------

clc;clear;close all;
%% Initialization
% The number of samples of all experiments
nSample = 1000;
% The number of realizations
nReal = 100;
% The delay
delay = 1;
% Leakage coefficient
leakage = 0;
% The MA process parameters
% MA process coefficents
MA_coeff = 0.9;
% The AR process order
order = length(MA_coeff);
% The noise power
varNoise = 0.5;

% GASS parameters
% The algorithm name
algo = {'Benveniste'};
% The initial time step
step0_GASS = 0.1;
% The learning rate for step size
rho_GASS = 5*1e-2;
% The GASS error
errorBen = cell(1,nReal);
% The GASS weight
weight_Ben = cell(1,nReal);

% The GNGD algorithm
% The initial time step
step0_GNGD = 1;
% The learning rate for step size
rho_GNGD = 5*1e-2;
% The GNGD error
errorGNGD = cell(1,nReal);
% The GNGD weight
weight_GNGD = cell(1,nReal);

%% Implement GASS and GNGD
for iReal = 1:nReal
   % Genrate MA process
   inputSig = random('Normal', 0, varNoise, nSample, 1).';
   desireSig = filter([1,0.9],1,inputSig);
   % The desired signal has a unit delay
   desireSig = [0,desireSig(1:end-1)];
   % GASS
  [weight_Ben{iReal}, errorBen{iReal}, ~] = funGASS(inputSig,desireSig,step0_GASS,order+1,delay,rho_GASS,leakage,string(algo{1}));
  % GNGD
  [weight_GNGD{iReal},errorGNGD{iReal},~] = funGNGD(inputSig,desireSig,order+1,step0_GNGD,rho_GNGD,delay,leakage);
end
av_errorBen = mean(cat(3,errorBen{:}).^2,3);
wError_Ben = MA_coeff - mean(cat(3, weight_Ben{:}),3); 
av_errorGNGD = mean(cat(3,errorGNGD{:}).^2,3);
wError_GNGD = MA_coeff - mean(cat(3, weight_GNGD{:}),3); 

%% Plot the result
figure;
subplot(1,2,1);
plot(wError_Ben(2,:),'LineWidth',2);
hold on;
plot(wError_GNGD(2,:),'LineWidth',2);
hold off;
xlabel('Number of iterations');
ylabel('Weight error');
title('Weight error curves for fixed and adaptive step size LMS');
set(gca,'fontsize',12);
grid on; grid minor;
legend('GASS','GNGD');
xlim([0,100]);
ylim([0,1]);

subplot(1,2,2);
plot(pow2db(abs(av_errorBen)),'LineWidth',2);
hold on
plot(pow2db(abs(av_errorGNGD)),'LineWidth',2);
hold off;
xlabel('Number of iterations');
ylabel('Averaged squared error (dB)');
title('The squared error curve for fixed and adaptive LMS');
set(gca,'fontsize',12);
grid on; grid minor;
legend('GASS','GNGD');
xlim([0,400]);