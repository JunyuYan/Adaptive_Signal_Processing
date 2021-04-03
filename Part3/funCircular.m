function [Cir_quotient, Cir_coeff] = funCircular(inputSig)
%% -------------------------------------------------------------------------
% This function is used to obtain the circularity quotient and circularity
% coefficent of a given signal 
% Input: inputSig- Input signal
% Output: Cir_quotient- Circularity quotient of the signal
%         Cir_coeff- Circularity coefficient of the signal

%% -------------------------------------------------------------------------
    % Covariance
    cov = mean(abs(inputSig).^2);
    % Pseudocovariance
    pseCov = mean(inputSig.^2);
    % Circularity quotient
    Cir_quotient = pseCov/cov;
    % Circularity coefficient
    Cir_coeff = abs(pseCov)/cov;
end