%% Task 2.1.a GASS algorithms
% Author: Junyu Yan
%--------------------------------------------------------------------------

clc;clear;close all;
%% Initialization
% The number of samples of all experiments
nSample = 1000;
% The number of realizations
nReal = 100;
% The step size of all experiments
step = [0.01 ; 0.1];
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

% LMS parameters
% The LMS error
errorLMS = cell(length(step),nReal);
% Average LMS error 
av_errorLMS = cell(length(step),1);
% The LMS weight
weight_LMS = cell(length(step),nReal);
% LMS weight error
wError_LMS = cell(length(step),1);

% GASS parameters
% The algorithm name
algo = {'Benveniste', 'Ang', 'Matthews'};
% The learning coefficient for Ang & Farhang
alpha = 0.8;
% The initial time step
step0 = 0.2;
% The learning rate for step size
rho = 5*1e-3;
% The GASS error
errorBen = cell(1,nReal);
errorAng = cell(1,nReal);
errorMatt = cell(1,nReal);
% The GASS weight
weight_Ben = cell(1,nReal);
weight_Ang = cell(1,nReal);
weight_Matt = cell(1,nReal);


%% LMS
for iStep = 1:length(step)
    for iReal = 1:nReal
        % Generate MA process
        inputSig = random('Normal', 0, varNoise, nSample, 1).';
        desireSig = filter([1,0.9],1,inputSig);
        % The desired signal has a unit delay
        desireSig = [0,desireSig(1:end-1)];
        % LMS
        [weight_LMS{iStep,iReal},errorLMS{iStep,iReal},~] = funLMS(inputSig,desireSig,order+1,step(iStep),delay,leakage);
    end
    av_errorLMS{iStep} = mean(cat(3,errorLMS{iStep,:}).^2,3);
    wError_LMS{iStep} = MA_coeff - mean(cat(3, weight_LMS{iStep,:}),3); 
end

%% GASS
for iReal = 1:nReal
    inputSig = random('Normal', 0, varNoise, nSample, 1).';
    desireSig = filter([1,0.9],1,inputSig);
   % The desired signal has a unit delay
   desireSig = [0,desireSig(1:end-1)];
   % GASS
   [weight_Ben{iReal}, errorBen{iReal}, ~] = funGASS(inputSig,desireSig,step0,order+1,delay,rho,leakage,string(algo{1}));
   [weight_Ang{iReal}, errorAng{iReal}, ~] = funGASS(inputSig,desireSig,step0,order+1,delay,rho,leakage,string(algo{2}),alpha);
   [weight_Matt{iReal}, errorMatt{iReal}, ~] = funGASS(inputSig,desireSig,step0,order+1,delay,rho,leakage,string(algo{3}));
end
av_errorBen = mean(cat(3,errorBen{:}).^2,3);
av_errorAng = mean(cat(3,errorAng{:}).^2,3);
av_errorMatt = mean(cat(3,errorMatt{:}).^2,3);
wError_Ben = MA_coeff - mean(cat(3, weight_Ben{:}),3); 
wError_Ang = MA_coeff - mean(cat(3, weight_Ang{:}),3); 
wError_Matt = MA_coeff - mean(cat(3, weight_Matt{:}),3); 

%% Plot the result
fig = figure;
subplot(1,2,1);
for iStep = 1:length(step)
    plot(wError_LMS{iStep}(2,:),'LineWidth',2);
    hold on;
end
plot(wError_Ben(2,:),'LineWidth',2);
hold on;
plot(wError_Ang(2,:),'LineWidth',2);
hold on;
plot(wError_Matt(2,:),'LineWidth',2);
hold off;
xlabel('Number of iterations');
ylabel('Weight error');
title('Weight error curves for fixed and adaptive step size LMS');
set(gca,'fontsize',12);
grid on; grid minor;
legend('\mu = 0.01','\mu = 0.1','Benveniste','Ang','Matthews','NumColumns',2);

subplot(1,2,2);
for iStep = 1:length(step)
    plot(pow2db(av_errorLMS{iStep}),'LineWidth',2);
    hold on;
end
plot(pow2db(av_errorBen),'LineWidth',2);
hold on
plot(pow2db(av_errorAng),'LineWidth',2);
hold on;
plot(pow2db(av_errorMatt),'LineWidth',2);
xlabel('Number of iterations');
ylabel('Averaged squared error (dB)');
title('The squared error curve for fixed and adaptive LMS');
set(gca,'fontsize',12);
grid on; grid minor;
legend('\mu = 0.01','\mu = 0.1','Benveniste','Ang','Matthews','NumColumns',2);