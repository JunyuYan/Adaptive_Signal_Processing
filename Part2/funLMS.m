function [weight_LMS, error_LMS, pred_LMS] = funLMS(inputSig,desireSig,order,step,delay,leakage)
%% -------------------------------------------------------------------------
% This function is used to implement the LMS algorithm
% Input: inputSig- The previous samples used for the prediction of future
%                  samples
%        desireSig- Desired output signal
%        step- Step size
%        order- The order of AR filter
%        delay- Delay in samples
%        leakage- leakage coefficient
% Output: weight_LMS-Weight of LMS
%         error_LMS-Error of LMS
%         pred_LMS-Predicted output by LMS
%% -------------------------------------------------------------------------
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
    % Define parameters
    % The number of input samples
    N = size(inputSig,2);
    % The weight
    weight = zeros(order,N+1);
    % The error
    error_LMS = zeros(1,N);
    % The predicted output
    pred_LMS = zeros(1,N);
    % The processed signal
    xn_LMS = zeros(order,N);
  
    % Obtain the processed signal x(n-m)
    for m = 1:order
        xn_LMS(m,:) = [zeros(1,m+delay-1),inputSig(1,1:(N-m-delay+1))];
    end
    
    % Iterate
    for i = 1:N
        pred_LMS(:,i) = weight(:,i).'* xn_LMS(:,i);
        error_LMS(:,i) = desireSig(:,i) - pred_LMS(:,i);
        weight(:,i+1) = (1-step*leakage)*weight(:,i) + step * error_LMS(:,i)*xn_LMS(:,i);
    end
    
    weight_LMS = weight(:,2:end);
end