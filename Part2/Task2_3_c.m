%% Task 2.3.c Compare performance between ANC and ALE
% Author: Junyu Yan
%--------------------------------------------------------------------------

clc;clear;close all;
%% Initialization
% The number of samples of all experiments
nSample = 1000;
% Normalized sampling frequency
fSample = 1;
% The number of realizations
nReal = 100;
% Learning step size
step = 0.01;
% Leakage
leakage = 0;
% Delays
delay = 3;
% Steady state offset
t0 = 50;
% The MA process parameters
% MA process coefficents
% MA process coefficents
MA_b = [1,0,0.5];
MA_a = 1;
% The noise power
varNoise = 1;
% Filter length
orderFilter = 6;

% The sinusoid parameters
% Time
t = (0:nSample-1)/fSample;
% Amplitude
aSine = 1;
% Normalized frequency
fSine = 5*1e-3;
% Sinusiod wave
xn = aSine*sin(2*pi*fSine*t);

% Paramters initialization
% Noise corrupted signal
sn = zeros(nReal,nSample);
% ALE Predicted output
pred_ALE = zeros(nReal,nSample);
% ALE MSE
SE_ALE = zeros(nReal,nSample-t0);
% ANC predicted noise
noise_ANC = zeros(nReal,nSample);
% ANC predicted signal
pred_ANC = zeros(nReal,nSample);
% ANC SE
SE_ANC = zeros(nReal,nSample-t0);

%% Prediced xn from ANC and ALE

for iReal = 1:nReal
    % Noise corrupted signal
    % Guassian noise
    vNoise = random('Normal', 0, varNoise, nSample, 1);
    Col_noise = filter(MA_b,MA_a,vNoise).';
    % Secondary noise
    noiseSec = 0.9*Col_noise+0.05;
    % Noise corrupted signal
    sn(iReal,:) = xn+Col_noise;
    inputSig_ALE = sn(iReal,:);
    desireSig_ALE = sn(iReal,:);
    % Apply ALE algorithm
    [~,~,pred_ALE(iReal,:)] = funLMS(inputSig_ALE,desireSig_ALE,orderFilter,step,delay,leakage);
    % ALE SE error
    SE_ALE(iReal,:) = (xn(t0+1:end)-pred_ALE(iReal,t0+1:end)).^2;
    
    inputSig_ANC = noiseSec;
    desireSig_ANC = [0,sn(iReal,1:end-1)];
    % Apply ANC algorithm
    [~,~,noise_ANC(iReal,:)] = funLMS(inputSig_ANC,desireSig_ANC,orderFilter,step,1,leakage);
    % Predicted xn
    pred_ANC(iReal,:) = desireSig_ANC - noise_ANC(iReal,:);
    % ANC SE error
    SE_ANC(iReal,:) = (xn(t0+1:end)-pred_ANC(iReal,t0+1:end)).^2;
end
% MSPE error
MSPE_ALE = mean(SE_ALE(:));
MSPE_ANC = mean(SE_ANC(:));

%% Plot results
% Plot the noise corrupted signal, the clean signal and predicted signal
figure;
subplot(1,2,1);
for iReal = 1:nReal
    % Plot noisy signal
    fig1 = plot(t,sn(iReal,:),'b','LineWidth',2);
    hold on;
end
for iReal = 1:nReal
    % Plot prediced signal
    fig2 = plot(t,pred_ALE(iReal,:),'r','LineWidth',2);
    hold on;
end
% Plot clean signal
fig3 = plot(t,xn,'k','LineWidth',2);
hold off;
title(sprintf('ALE signal, MSPE = %.2f dB', pow2db(MSPE_ALE)));
xlabel('Time (Samples)');
ylabel('Amplitude');
legend([fig1,fig2,fig3],'Noisy','ALE','Clean','NumColumns',3);
ylim([-5,5]);
grid on; grid minor;

subplot(1,2,2);
for iReal = 1:nReal
    % Plot noisy signal
    fig1 = plot(t,sn(iReal,:),'b','LineWidth',2);
    hold on;
end
for iReal = 1:nReal
    % Plot prediced signal
    fig2 = plot(t,pred_ANC(iReal,:),'r','LineWidth',2);
    hold on;
end
% Plot clean signal
fig3 = plot(t,xn,'k','LineWidth',2);
hold off;
title(sprintf('ANC signal, MSPE = %.2f dB', pow2db(MSPE_ANC)));
xlabel('Time (Samples)');
ylabel('Amplitude');
legend([fig1,fig2,fig3],'Noisy','ANC','Clean','NumColumns',3);
ylim([-5,5]);
grid on; grid minor;

% Plot the mean 100 realizations signal for ANC and ALE
figure;
plot(t,mean(pred_ALE),'b','LineWidth',2);
hold on;
plot(t,mean(pred_ANC),'r','LineWidth',2);
hold on;
plot(t,xn,'k','LineWidth',2);
title(sprintf('ALE vs ANC: Ensemble mean'));
xlabel('Time (Samples)');
ylabel('Amplitude');
legend('ALE','ANC','Clean','NumColumns',3);
grid on; grid minor;