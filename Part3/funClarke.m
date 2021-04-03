function [Clarke_V] = funClarke(inputV)
%% -------------------------------------------------------------------------
% This function is used to apply clarke transform to a given signal
% Input: inputV- Input Voltage
% Output: Clarke_V- Clarke voltage

%% -------------------------------------------------------------------------
    % Clarke matrix
    C = sqrt(2/3)*[sqrt(2)/2 sqrt(2)/2 sqrt(2)/2; 1 -1/2 -1/2; 0 sqrt(3)/2 -sqrt(3)/2];
    % Clarke transform
    Clarke_tmp = C * inputV;
    % Clarke voltage
    Clarke_V = complex(Clarke_tmp(2,:),Clarke_tmp(3,:));
end