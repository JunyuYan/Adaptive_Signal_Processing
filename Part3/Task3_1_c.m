%% Task 3.1.c Circularity of balanced and unbalanced system
% Author: Junyu Yan
%--------------------------------------------------------------------------

clc;clear;close all;
%% Initialization
% The number of samples
nSample = 1000;
% The number of phase shift
nPhase = 3;
% The sampling frequency
fs = 1e4;
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
%% Balanced three phase power systems
% Balanced three phase
balanced_V = Amp .* cos(2*pi*fo/fs*t+delta+phase);
% Clake voltage
Clarke_bal = funClarke(balanced_V);
% Circularity
[~,Cir_bal] = funCircular(Clarke_bal);

%% Unbalanced three phase power system
% Unbalanced amplitude
% Amplitude difference
ampDiff = 0.2:0.2:0.8;
% Clarke voltage
Clarke_unAmp = zeros(length(ampDiff),nSample);
% Circularity
Cir_unAmp = zeros(length(ampDiff),1);
for iAmp = 1:length(ampDiff)
    % Unbalanced three phase system
    unAmp_V = (Amp+[-ampDiff(iAmp);0;ampDiff(iAmp)]).*cos(2*pi*fo/fs*t+delta+phase);
    % Clake voltage
    Clarke_unAmp(iAmp,:) = funClarke(unAmp_V);
    % Circularity
    [~,Cir_unAmp(iAmp)] = funCircular(Clarke_unAmp(iAmp,:));
end

% Unbalanced phases
% Phases difference
phiDiff = 0.05:0.05:0.2;
% Clarke voltage
Clarke_unPhi = zeros(length(phiDiff),nSample);
% Circularity
Cir_unPhi = zeros(length(phiDiff),1);
for iPhi = 1:length(phiDiff)
    % Phase shift
    delta = zeros(3,1)+[0;-phiDiff(iPhi)*pi;phiDiff(iPhi)*pi];
    % Unbalanced three phase system
    unPhi_V = Amp.*cos(2*pi*fo/fs*t+delta+phase);
    % Clake voltage
    Clarke_unPhi(iPhi,:) = funClarke(unPhi_V);
    % Circularity
    [~,Cir_unPhi(iPhi)] = funCircular(Clarke_unPhi(iPhi,:));
end
%% Result plot
figure;
Color = get(groot, 'factoryAxesColorOrder');
% Balanced three phase
scatter(real(Clarke_bal),imag(Clarke_bal),30,Color(1,:),'filled');
title(sprintf('Balanced circularity diagram with |\\rho| = %.2f', Cir_bal));
grid on; grid minor;
xlabel('Real part');
ylabel('Imaginary part');
set(gca,'fontsize',12);
axis square;
figure;
% Unbalanced three phase with diferent amplitudes
legendstr = cell(length(ampDiff),1);
subplot(1,2,1);
for iAmp = 1:length(ampDiff)
    scatter(real(Clarke_unAmp(iAmp,:)),imag(Clarke_unAmp(iAmp,:)),30,Color(iAmp+1,:),'filled');
    legendstr{iAmp} = sprintf('\\DeltaV = %.1f, \\rho = %.2f', ampDiff(iAmp),Cir_unAmp(iAmp));
    hold on;
    title(sprintf('Circularity diagram with unbalanced amplitudes'));
    grid on; grid minor;
    xlabel('Real part');
    ylabel('Imaginary part');
    set(gca,'fontsize',12);
end
legend(legendstr);
xlim([-2,2]);
ylim([-2,2]);
axis square;
% Unbalanced three phase with diferent phases
legendstr = cell(length(phiDiff),1);
subplot(1,2,2);
for iPhi = 1:length(phiDiff)
    scatter(real(Clarke_unPhi(iPhi,:)),imag(Clarke_unPhi(iPhi,:)));
    legendstr{iPhi} = sprintf('\\Delta\\phi = %.2f\\pi, \\rho = %.2f', phiDiff(iPhi),Cir_unPhi(iPhi));
    hold on;
    title(sprintf('Circularity diagram with unbalanced phases'));
    grid on; grid minor;
    xlabel('Real part');
    ylabel('Imaginary part');
    set(gca,'fontsize',12);
end
legend(legendstr);
xlim([-2,2]);
ylim([-2,2]);
axis square;