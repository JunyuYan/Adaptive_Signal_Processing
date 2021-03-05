%% Task 1.3.e Correlation Estimation
% This task is to find the desired line spectra using MUSIC method
% Author: Junyu Yan
%--------------------------------------------------------------------------


clc;clear;close all;
%% Initialization
% The number of samples
N = [20,30,40,50];
% The normalized sampling frequency
fs = 1;
% The frequency of sinusoidal signal
f_sin = [0.3, 0.32];
% The amplitude of sinusoidal signal
A_sin = [1, 1];
% Noise power
nPower = 0.2;
% The number of dft
nDFT = 512;
% The figure index
Fig_index = 0;
%The number of realisation
nRe = 100;
% The PSD matrix is 
psd_matrix = zeros(nDFT,nRe,length(N));
%% Generate and plot the PSD for 100 realisations
for i = 1:length(N)
    Fig_index = Fig_index+1;
    figure(1);
    subplot(length(N),1,Fig_index);
    for iRe = 1:nRe
        % The sampling time
        n = (0: N(i) - 1) / fs;
        % Generate noise
        noise = sqrt(nPower/2)*(randn(size(n))+1j*randn(size(n)));
        % Generate the signal
        x = A_sin(1)*exp(1j*2*pi*f_sin(1)*n)+A_sin(2)*exp(1j*2*pi*f_sin(2)*n)+noise;
        % Music estimate
        [X,R] = corrmtx(x,14,'mod');
        [S,F] = pmusic(R,2,nDFT,fs,'corr');
        % Plot the result
        F1 = plot(F,S,'c','linewidth',1);
        hold on;
        psd_matrix(:,iRe,i) = S; 
    end
    %% Calculate and plot the mean PSD
    psd_mean = mean(psd_matrix(:,:,i),2);
    F2 = plot(F,psd_mean,'b','linewidth',2);
    set(gca,'xlim',[0.25, 0.40]);
    grid on;grid minor;
    xlabel('Hz');
    ylabel('Pseudospectrum');
    title({['PSD with MUSIC method N=',num2str(N(i))]});
    legend([F1,F2],'Realisation','Mean');
end
%% Generate and plot the standard derivation
figure(2);
Fig_index = 0;
for i = 1:length(N)
    Fig_index = Fig_index+1;
    subplot(length(N),1,Fig_index);
    psd_std = std(psd_matrix(:,:,i),[],2);
    plot(F,psd_std,'r','linewidth',2);
    set(gca,'xlim',[0.25, 0.40]);
    grid on;grid minor; 
    xlabel('Hz');
    ylabel('Std');
    title({['Standard Derivative with MUSIC method N=',num2str(N(i))]});
end


