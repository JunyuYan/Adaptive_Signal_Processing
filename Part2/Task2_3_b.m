%% Task 2.3.a Order vs MSPE for adaptive line enhancer
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
delay = 1:25;
% Steady state offset
t0 = 50;
% Filter length
orderFilter = 1:20;
% The MA process parameters
% MA process coefficents
MA_b = [1,0,0.5];
MA_a = 1;
% The noise power
varNoise = 1;

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
sn = cell(length(delay),nReal);
% Predicted output
pred_ALE = cell(length(delay),nReal);
% MSE
MSE_ALE = cell(length(delay),nReal);
% MSPE
MSPE_ALE = zeros(length(delay),length(orderFilter));


%% Calculate MSPE with different delay and filter order for ALE
for iOrder = orderFilter
    for iDelay = delay
        for iReal = 1:nReal
            % Noise corrupted signal
            % Guassian noise
            vNoise = random('Normal', 0, varNoise, nSample, 1);
            Col_noise = filter(MA_b,MA_a,vNoise).';
            sn{iDelay,iReal} = xn+Col_noise;
            inputSig = sn{iDelay,iReal};
            desireSig = sn{iDelay,iReal};
            % Apply LMS algorithm
            [~,~,pred_ALE{iDelay,iReal}] = funLMS(inputSig,desireSig,iOrder,step,iDelay,leakage);
            % MSE error
            MSE_ALE{iDelay,iReal} = (xn(t0+1:end)-pred_ALE{iDelay,iReal}(t0+1:end)).^2;
        end
        % MSPE error
        MSPE_ALE(iDelay,iOrder) = mean(cell2mat(MSE_ALE(iDelay,:)));
    end
end

%% Plot results
% Plot the MSPE against filter orders and delay
figure; 
subplot(1,3,1);
for iOrder = 5:5:20
    plot(pow2db(MSPE_ALE(:,iOrder)),'LineWidth',2);
    hold on;
end
xlabel('Delay (Sample)');
ylabel('MSPE (dB)');
xlim([1,25]);
xticks(1:1:25);
title(sprintf('ALE:MSPE against filter orders and delay'));
legend('M = 5','M = 10','M = 15','M = 20');
grid on; grid minor;
% Plot the MSPE against delay
subplot(1,3,2);
plot(pow2db(MSPE_ALE(:,5)),'LineWidth',2);
xlabel('Delay (Sample)');
ylabel('MSPE (dB)');
xlim([1,25]);
xticks(1:1:25);
title(sprintf('ALE:MSPE against \\Delta, M = %d',5));
grid on; grid minor;
% Plot the MSPE against filter order
subplot(1,3,3);
plot(pow2db(MSPE_ALE(3,:).'),'LineWidth',2);
xlabel('Filter Order');
ylabel('MSPE (dB)');
title(sprintf('ALE:MSPE against M, \\Delta = %d',3));
xlim([1,20]);
xticks(1:1:20);
grid on; grid minor;

