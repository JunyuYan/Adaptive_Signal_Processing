function [weight_DFT, error_DFT, pred_DFT] = funDFT_CLMS(inputSig,desireSig,step,leakage)
%% -------------------------------------------------------------------------
% This function is used to implement the DFT_CLMS algorithm
% Input: inputSig- The previous samples used for the prediction of future
%                  samples
%        desireSig- Desired output signal
%        step- Step size
%        leakage- leakage coefficient
% Output: weight_DFT-Weight of DFT-CLMS
%         error_DFT-Error of DFT-CLMS
%         pred_DFT-Predicted output by DFT-CLMS
%% -------------------------------------------------------------------------
    % Check the input
    if ~ismatrix(inputSig)
        error('The input signal should be a matrix');
    end
    if ~isvector(desireSig)
        error('The desired signal should be a vector');
    end
    if ~isscalar(step)
        error('The step size should be a scalar');
    end
    if ~isscalar(leakage)
        error('The leakage should be a scalar');
    end
    % Define parameters
    % The number of input samples
    [M,N] = size(inputSig);
    % The weight
    weight = zeros(M,N+1);
    % The error
    error_DFT = zeros(1,N);
    % The predicted output
    pred_DFT = zeros(1,N);
    % The processed signal
    xnDFT_CLMS = inputSig;
    
    % Iterate
    for i = 1:N
        pred_DFT(:,i) = weight(:,i)'* xnDFT_CLMS(:,i);
        error_DFT(:,i) = desireSig(:,i) - pred_DFT(:,i);
        weight(:,i+1) = (1-step*leakage)*weight(:,i) + step * conj(error_DFT(:,i))*xnDFT_CLMS(:,i);
    end
    
    weight_DFT = weight(:,2:end);
end