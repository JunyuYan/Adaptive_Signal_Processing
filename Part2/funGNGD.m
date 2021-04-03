function [weight_GNGD, error_GNGD, pred_GNGD] = funGNGD(inputSig,desireSig,order,step0,rho,delay,leakage)
%% -------------------------------------------------------------------------
% This function is used to implement the LMS algorithm
% Input: inputSig- The previous samples used for the prediction of future
%                  samples
%        desireSig- Desired output signal
%        step0- Initial step size
%        order- The order of AR filter
%        rho- The learning rate for regularization factor
%        delay- Delay in samples
%        leakage- leakage coefficient
% Output: weight_GNGD-Weight of GNGD
%         error_GNGD-Error of GNGD
%         pred_GNGD-Predicted output by GNGD
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
    if ~isscalar(step0)
        error('The step size should be a scalar');
    end
    if ~isscalar(rho)
        error('The learning rate should be a scalar');
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
    error_GNGD = zeros(1,N);
    % The predicted output
    pred_GNGD = zeros(1,N);
    % The processed signal
    xn_GNGD = zeros(order,N);
    % The regularization factor
    reg_factor = ones(1,N+1)./step0;
  
    % Obtain the processed signal x(n-m)
    for m = 1:order
        xn_GNGD(m,:) = [zeros(1,m+delay-1),inputSig(1,1:(N-m-delay+1))];
    end
    
    % Iterate
    for n = 1:N
        pred_GNGD(n) = weight(:,n).'* xn_GNGD(:,n);
        error_GNGD(n) = desireSig(:,n) - pred_GNGD(:,n);
        weight(:,n+1) = (1-step0*leakage)*weight(:,n) + 1/(reg_factor(n)+xn_GNGD(:,n).'*xn_GNGD(:,n))*error_GNGD(:,n)*xn_GNGD(:,n);
        if n>1
            reg_factor(n+1) = reg_factor(n)-rho*step0*(error_GNGD(n)*error_GNGD(n-1)*xn_GNGD(:,n).'*xn_GNGD(:,n-1))/(reg_factor(n-1)+xn_GNGD(:,n-1)'*xn_GNGD(:,n-1))^2;
        end
    end
    weight_GNGD = weight(:,2:end);
end