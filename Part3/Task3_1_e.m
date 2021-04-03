%% Task 3.1.e CLMS vs ACLMS for three phase power systems
% Author: Junyu Yan
%--------------------------------------------------------------------------

clc;clear;close all;
%% Initialization
% The number of samples
nSample = 1000;
% The number of phase shift
nPhase = 3;
% The sampling frequency
fs = 1000;
% The system frequency
fo = 50;
% Time
t = 0:nSample-1;
% Initial phase
phase = [0;-2*pi/3;2*pi/3];
% The amplitude 
Amp = ones(nPhase,1);
% Phase shift
delta = zeros(nPhase,1);
% CLMS and ACLMS parameters
% Learning step size
step = 0.05;
% Leakage
leakage = 0;
% Delays
delay = 1;
% The order
orderFilter = 1;
% Paramters initialization

%% Balanced three phase power systems
% Balanced three phase
balanced_V = Amp .* cos(2*pi*fo/fs*t+delta+phase);
% Clake voltage
Clarke_bal = funClarke(balanced_V);
% Circularity
[~,Cir_bal] = funCircular(Clarke_bal);
% CLMS and ACLMS algorithm
inputSig = Clarke_bal;
desireSig = Clarke_bal;
% Apply CLMS algorithm
[wBal_CLMS,eBal_CLMS,~] = funCLMS(inputSig,desireSig,orderFilter,step,delay,leakage);
% Apply CLMS algorithm
[hBal_ACLMS,gBal_ACLMS,eBal_ACLMS,~] = funACLMS(inputSig,desireSig,orderFilter,step,delay,leakage);
% Frequency
foBal_CLMS = fs/(2*pi)*atan(imag(wBal_CLMS)./real(wBal_CLMS));
foBal_ACLMS = fs/(2*pi)*atan(sqrt(imag(hBal_ACLMS).^2-abs(gBal_ACLMS).^2)./real(hBal_ACLMS));
% Plot the learning curve
figure;
subplot(1,2,1);
plot(pow2db(abs(eBal_CLMS).^2),'LineWidth',2);
hold on;
plot(pow2db(abs(eBal_ACLMS).^2),'LineWidth',2);
grid on;grid minor;
title('Learning curve for balanced system');
xlabel('Time (Sample)');
ylabel('Squared error (dB)');
legend('CLMS','ACLMS');
set(gca,'fontsize',10);
% Plot the frequency estimation
subplot(1,2,2);
plot(abs(foBal_CLMS),'LineWidth',2);
hold on;
plot(abs(foBal_ACLMS),'LineWidth',2);
hold on;
plot([0,nSample],[50,50],'--','color','k','LineWidth',2);
grid on;grid minor;
title('Nominal frequency estimation for balanced system');
xlabel('Time (Sample)');
ylabel('Frequency estimation (Hz)');
legend('CLMS','ACLMS','True');
set(gca,'fontsize',10);
ylim([0,100]);

%% Unbalanced amplitude three phase power systems
% Amplitude difference
ampDiff = 0.6;
% Phase difference
phiDiff = 0.1;
% Phase shift
delta = zeros(3,1)+[0;-phiDiff*pi;phiDiff*pi];
% Unbalanced three phase
unAmp_V = (Amp+[-ampDiff;0;ampDiff]).* cos(2*pi*fo/fs*t+delta+phase);
% Clake voltage
Clarke_unAmp = funClarke(unAmp_V);
% Circularity
[~,Cir_unAmp] = funCircular(Clarke_unAmp);
% CLMS and ACLMS algorithm
inputSig = Clarke_unAmp;
desireSig = Clarke_unAmp;
% Apply CLMS algorithm
[wUnamp_CLMS,eUnamp_CLMS,~] = funCLMS(inputSig,desireSig,orderFilter,step,delay,leakage);
% Apply CLMS algorithm
[hUnamp_ACLMS,gUnamp_ACLMS,eUnamp_ACLMS,~] = funACLMS(inputSig,desireSig,orderFilter,step,delay,leakage);
% Frequency
foUnamp_CLMS = fs/(2*pi)*atan(imag(wUnamp_CLMS)./real(wUnamp_CLMS));
foUnamp_ACLMS = fs/(2*pi)*atan(sqrt(imag(hUnamp_ACLMS).^2-abs(gUnamp_ACLMS).^2)./real(hUnamp_ACLMS));
% Plot the learning curve
figure;
subplot(1,2,1);
plot(pow2db(abs(eUnamp_CLMS).^2),'LineWidth',2);
hold on;
plot(pow2db(abs(eUnamp_ACLMS).^2),'LineWidth',2);
grid on;grid minor;
title('Learning curve for unbalanced amplitude system');
xlabel('Time (Sample)');
ylabel('Squared error (dB)');
legend('CLMS','ACLMS');
set(gca,'fontsize',10);
% Plot the frequency estimation
subplot(1,2,2);
plot(abs(foUnamp_CLMS),'LineWidth',2);
hold on;
plot(abs(foUnamp_ACLMS),'LineWidth',2);
hold on;
plot([0,nSample],[50,50],'--','color','k','LineWidth',2);
grid on;grid minor;
title('Nominal frequency estimation for unbalanced amplitude system');
xlabel('Time (Sample)');
ylabel('Frequency estimation (Hz)');
legend('CLMS','ACLMS','True');
set(gca,'fontsize',10);
ylim([0,100]);

%% Unbalanced amplitude three phase power systems
% Phase difference
phiDiff = 0.1;
% Phase shift
delta = zeros(3,1)+[0;-phiDiff*pi;phiDiff*pi];
% Unbalanced three phase
unPhi_V = Amp.* cos(2*pi*fo/fs*t+delta+phase);
% Clake voltage
Clarke_unPhi = funClarke(unPhi_V);
% Circularity
[~,Cir_unPhi] = funCircular(Clarke_unPhi);
% CLMS and ACLMS algorithm
inputSig = Clarke_unPhi;
desireSig = Clarke_unPhi;
% Apply CLMS algorithm
[wUnphi_CLMS,eUnphi_CLMS,~] = funCLMS(inputSig,desireSig,orderFilter,step,delay,leakage);
% Apply CLMS algorithm
[hUnphi_ACLMS,gUnphi_ACLMS,eUnphi_ACLMS,~] = funACLMS(inputSig,desireSig,orderFilter,step,delay,leakage);
% Frequency
foUnphi_CLMS = fs/(2*pi)*atan(imag(wUnphi_CLMS)./real(wUnphi_CLMS));
foUnphi_ACLMS = fs/(2*pi)*atan(sqrt(imag(hUnphi_ACLMS).^2-abs(gUnphi_ACLMS).^2)./real(hUnphi_ACLMS));
% Plot the learning curve
figure;
subplot(1,2,1);
plot(pow2db(abs(eUnphi_CLMS).^2),'LineWidth',2);
hold on;
plot(pow2db(abs(eUnphi_ACLMS).^2),'LineWidth',2);
grid on;grid minor;
title('Learning curve for unbalanced phase system');
xlabel('Time (Sample)');
ylabel('Squared error (dB)');
legend('CLMS','ACLMS');
set(gca,'fontsize',10);
% Plot the frequency estimation
subplot(1,2,2);
plot(abs(foUnphi_CLMS),'LineWidth',2);
hold on;
plot(abs(foUnphi_ACLMS),'LineWidth',2);
hold on;
plot([0,nSample],[50,50],'--','color','k','LineWidth',2);
grid on;grid minor;
title('Nominal frequency estimation for unbalanced phase system');
xlabel('Time (Sample)');
ylabel('Frequency estimation (Hz)');
legend('CLMS','ACLMS','True');
set(gca,'fontsize',10);
ylim([0,100]);