function [h_ACLMS, g_ACLMS, error_ACLMS, pred_ACLMS] = funACLMS(inputSig,desireSig,order,step,delay,leakage)
%% -------------------------------------------------------------------------
% This function is used to implement the ACLMS algorithm
% Input: inputSig- The previous samples used for the prediction of future
%                  samples
%        desireSig- Desired output signal
%        step- Step size
%        order- The order of filter
%        delay- Delay in samples
%        leakage- leakage coefficient
% Output: h_ACLMS- h weight of ACLMS
%         g_ACLMS- g weight of ACLMS
%         error_ACLMS-Error of ACLMS
%         pred_ACLMS-Predicted output by ACLMS
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
    hWeight = zeros(order,N+1);
    gWeight = zeros(order,N+1);
    % The error
    error_ACLMS = zeros(1,N);
    % The predicted output
    pred_ACLMS = zeros(1,N);
    % The processed signal
    xn_ACLMS = zeros(order,N);
  
    % Obtain the processed signal x(n-m)
    for m = 1:order
        xn_ACLMS(m,:) = [zeros(1,m+delay-1),inputSig(1,1:(N-m-delay+1))];
    end
    
    % Iterate
    for n = 1:N
        pred_ACLMS(:,n) = hWeight(:,n)'* xn_ACLMS(:,n)+gWeight(:,n)'*conj(xn_ACLMS(:,n));
        error_ACLMS(:,n) = desireSig(:,n) - pred_ACLMS(:,n);
        hWeight(:,n+1) = (1-step*leakage)*hWeight(:,n) + step * conj(error_ACLMS(:,n))*xn_ACLMS(:,n);
        gWeight(:,n+1) = (1-step*leakage)*gWeight(:,n) + step * conj(error_ACLMS(:,n))*conj(xn_ACLMS(:,n));
    end
    h_ACLMS = hWeight(:,2:end);
    g_ACLMS = gWeight(:,2:end);
end