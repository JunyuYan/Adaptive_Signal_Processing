%% Task 3.1.a CLMS vs ALMS for WLMA model
% Author: Junyu Yan
%--------------------------------------------------------------------------

clc;clear;close all;
%% Initialization
% The number of samples of all experiments
nSample = 1000;
% The number of realizations
nReal = 100;
% CLMS and ACLMS parameters
% Learning step size
step = 0.1;
% Leakage
leakage = 0;
% Delays
delay = 1;

% The MA process parameters
% MA process coefficents
MA_coeff = [1.5+1i,2.5-0.5i];
% The noise power
varNoise = 1;
% The order
orderFilter = length(MA_coeff);

% Paramters initialization
% CLMS error
error_CLMS = zeros(nReal,nSample);
% ACLMS error
error_ACLMS = zeros(nReal,nSample);

%% Generate the WLMA process
% Generate the white noise
noise = sqrt(varNoise/2) * (randn(nReal,nSample)+1i*randn(nReal,nSample));
% Circular coefficient of white noise
[~,Cir_noise] = funCircular(noise);
% The input signal
yn = noise+MA_coeff(1)*[zeros(nReal,1),noise(:,1:end-1)]+MA_coeff(2)*[zeros(nReal,1),conj(noise(:,1:end-1))];
% Circular coefficient of WLMA signal
[Cir_yn,~] = funCircular(yn);
%% CLMS and ACLMS algorithm
for iReal = 1:nReal
    inputSig = noise(iReal,:);
    desireSig = [0,yn(iReal,1:end-1)];
    % Apply CLMS algorithm
    [~,error_CLMS(iReal,:),~] = funCLMS(inputSig,desireSig,orderFilter,step,delay,leakage);
    % Apply CLMS algorithm
    [~,~,error_ACLMS(iReal,:),~] = funACLMS(inputSig,desireSig,orderFilter,step,delay,leakage);
end
av_errorCLMS = mean(abs(error_CLMS).^2);
av_errorACLMS = mean(abs(error_ACLMS).^2);

%% Plot results
% Plot the white noise and WLMA signal
figure;
subplot(2,2,1);
scatter(real(noise(:)),imag(noise(:)),1,'r');
title(sprintf('White noise with |\\rho| = %.2f', Cir_noise));
grid on; grid minor;
xlabel('Real part');
ylabel('Imaginary part');
set(gca,'fontsize',10);
subplot(2,2,2);
scatter(real(yn(:)),imag(yn(:)),1,'b');
title(sprintf('WLMA signal with |\\rho| = %.2f', Cir_yn));
grid on; grid minor;
xlabel('Real part');
ylabel('Imaginary part');
set(gca,'fontsize',10);
subplot(2,1,2);
plot(pow2db(av_errorCLMS),'LineWidth',2);
hold on;
plot(pow2db(av_errorACLMS),'LineWidth',2);
title('Learning curves for CLMS and ACLMS');
xlabel('Time (Sample)');
ylabel('Squared error (dB)');
grid on; grid minor;
legend('CLMS','ACLMS');
set(gca,'fontsize',10);
