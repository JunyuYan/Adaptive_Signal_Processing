function [weight_CLMS, error_CLMS, pred_CLMS] = funCLMS(inputSig,desireSig,order,step,delay,leakage)
%% -------------------------------------------------------------------------
% This function is used to implement the CLMS algorithm
% Input: inputSig- The previous samples used for the prediction of future
%                  samples
%        desireSig- Desired output signal
%        step- Step size
%        order- The order of filter
%        delay- Delay in samples
%        leakage- leakage coefficient
% Output: weight_CLMS-Weight of CLMS
%         error_CLMS-Error of CLMS
%         pred_CLMS-Predicted output by CLMS
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
    error_CLMS = zeros(1,N);
    % The predicted output
    pred_CLMS = zeros(1,N);
    % The processed signal
    xn_CLMS = zeros(order,N);
  
    % Obtain the processed signal x(n-m)
    for m = 1:order
        xn_CLMS(m,:) = [zeros(1,m+delay-1),inputSig(1,1:(N-m-delay+1))];
    end
    
    % Iterate
    for i = 1:N
        pred_CLMS(:,i) = weight(:,i)'* xn_CLMS(:,i);
        error_CLMS(:,i) = desireSig(:,i) - pred_CLMS(:,i);
        weight(:,i+1) = (1-step*leakage)*weight(:,i) + step * conj(error_CLMS(:,i))*xn_CLMS(:,i);
    end
    
    weight_CLMS = weight(:,2:end);
end