function [weight_PM, error_PM, pred_PM] = funPerception(inputSig,desireSig,order,step,delay,leakage,scale,bias,weightInit)
%% -------------------------------------------------------------------------
% This function is used to implement the Perception model
% Input: inputSig- The previous samples used for the prediction of future
%                  samples
%        desireSig- Desired output signal
%        step- Step size
%        order- The order of AR filter
%        delay- Delay in samples
%        leakage- leakage coefficient
%        scale- Scaling factor for activation function
%        bias- Indicating whether there is bias
%        weightInit- Initial weight for training
% Output: weight_PM-Weight of perception model
%         error_PM-Error of perception model
%         pred_PM-Predicted output by perception model
%% -------------------------------------------------------------------------
    % Check the number of parameters
    if nargin < 9
        weightInit = [];
    elseif nargin < 8
        bias = 0;
    end
     
    % Check the input
    if ~isvector(inputSig)
        error('The input signal should be a vector');
    end
    if ~isvector(desireSig)
        error('The desired signal should be a vector');
    end
    if ~isscalar(order)
        error('The order should be a scalar');
    end
    if ~isscalar(step)
        error('The step size should be a scalar');
    end
    if ~isscalar(delay)
        error('The delay should be a scalar');
    end
    if ~isscalar(leakage)
        error('The leakage should be a scalar');
    end
    if ~isscalar(scale)
        error('The scale for activation function should be a scalar');
    end
    if ~isscalar(bias)
        error('The bias for activation function should be a scalar');
    end

    % Define parameters
    % The number of input samples
    N = size(inputSig,2);
    % The processed signal
    xn_LMS = zeros(order,N);
    % The error
    error_PM = zeros(1,N);
    % The predicted output
    pred_PM = zeros(1,N);
    
    % Obtain the processed signal x(n-m)
    for m = 1:order
        xn_LMS(m,:) = [zeros(1,m+delay-1),inputSig(1,1:(N-m-delay+1))];
    end
    % Bias signal
    if bias
        xn_LMS = [ones(1,size(xn_LMS,2)); xn_LMS];
    end
    
    % Define weight
    if isempty(weightInit)
        % New filter order
        M = size(xn_LMS,1);
        % The weight
        weight = zeros(M,N+1);
    else
        weight = weightInit;
    end
    
    % Iterate
    for i = 1:N
        pred_PM(:,i) = scale*tanh(weight(:,i).'* xn_LMS(:,i));
        error_PM(:,i) = desireSig(:,i) - pred_PM(:,i);
        weight(:,i+1) = (1-step*leakage)*weight(:,i) + (step * scale*(1-tanh(weight(:,i).'* xn_LMS(:,i))^2)*error_PM(:,i)).*xn_LMS(:,i);
    end
    
    weight_PM = weight(:,2:end);
end