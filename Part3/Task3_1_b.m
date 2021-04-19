%% Task 3.1.b CLMS and ACLMS for wind data
% Author: Junyu Yan
%--------------------------------------------------------------------------

clc;clear;close all;
load('../Data/wind-dataset/low-wind.mat');
wind(1,:) = complex(v_east,v_north);
load('../Data/wind-dataset/medium-wind.mat');
wind(2,:) = complex(v_east,v_north);
load('../Data/wind-dataset/high-wind.mat');
wind(3,:) = complex(v_east,v_north);
%% Initialization
% The number of samples
nSample = length(wind);
% The number of wind data
nWind = size(wind,1);
% CLMS and ACLMS parameters
% Learning step size
step = [0.1,0.005,0.001];
% Leakage
leakage = 0;
% Delays
delay = 1;
% The order
orderFilter = 1:1:25;

% Paramters initialization
% Circular coefficient
Cir_coeff = zeros(nWind,1);
% CLMS error
error_CLMS = cell(nWind,length(orderFilter));
% ACLMS error
error_ACLMS = cell(nWind,length(orderFilter));
% CLMS MSPE error
mspe_CLMS = zeros(nWind,length(orderFilter));
% ACLMS MSPE error
mspe_ACLMS = zeros(nWind,length(orderFilter));

%% Obtain the circularity of data
% Circular coefficient of wind data
for iWind = 1:nWind
    [~,Cir_coeff(iWind)] = funCircular(wind(iWind,:));
end

%% CLMS and ACLMS algorithm
for iWind = 1:nWind
    for iOrder = 1:length(orderFilter)
        inputSig = wind(iWind,:);
        desireSig = wind(iWind,:);
        % Apply CLMS algorithm
        [~,error_CLMS{iWind,iOrder},~] = funCLMS(inputSig,desireSig,iOrder,step(iWind),delay,leakage);
        % Apply ACLMS algorithm
        [~,~,error_ACLMS{iWind,iOrder},~] = funACLMS(inputSig,desireSig,iOrder,step(iWind),delay,leakage);
        % MSPE error
        mspe_CLMS(iWind,iOrder) = mean(abs(error_CLMS{iWind,iOrder}).^2);
        mspe_ACLMS(iWind,iOrder) = mean(abs(error_ACLMS{iWind,iOrder}).^2);
    end
end

%% Plot results
% Circularity plot and learning curve of wind data
label = {'Low','Medium','High'};
COLORS = get(groot, 'factoryAxesColorOrder');
for iWind = 1:nWind
    figure;
    subplot(1,2,1);
    scatter(real(wind(iWind,:)),imag(wind(iWind,:)),30,COLORS(iWind,:),'filled');
    title(sprintf('%s wind with |\\rho| = %.2f', label{iWind},Cir_coeff(iWind,:)));
    grid on; grid minor;
    xlabel('v_{east}');
    ylabel('v_{north}');
    set(gca,'fontsize',10);
    subplot(1,2,2);
    plot(pow2db(mspe_CLMS(iWind,:)),'LineWidth',2);
    hold on;
    plot(pow2db(mspe_ACLMS(iWind,:)),'LineWidth',2);
    title(sprintf('%s Wind: Learning curve for CLMS and ACLMS',label{iWind}));
    grid on; grid minor;
    xlabel('Model order');
    ylabel('MSPE (dB)');
    legend('CLMS','ACLMS');
    set(gca,'fontsize',10);
end

