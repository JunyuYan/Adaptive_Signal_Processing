function [weight_GASS, error_GASS, mu_GASS] = funGASS(inputSig,desireSig,step0,order,rho,leakage,algo,alpha)
%% -------------------------------------------------------------------------
% This function is used to implement the LMS algorithm
% Input: inputSig- The previous samples used for the prediction of future
%                  samples
%        desireSig- Desired output signal
%        step0- Initial step size
%        order- The order of AR process
%        rho- Learning rate fot the step size
%        leakage- leakage coefficient
%        algo- The algorithm name as {'Benveniste','Ang','Matthews'}
%        alpha- The adaptive coefficient for Ang & Farhang
% Output: weight_GASS- Weight of GASS
%         error_GASS- Error of GASS
%         mu_GASS- Predicted step size by GASS
%% -------------------------------------------------------------------------
    % Check the input
    if ~isvector(inputSig)
        error('The input signal should be a vector');
    end
    if ~isvector(desireSig)
        error('The desired signal should be a vector');
    end
    if ~isscalar(step0)
        error('The intial step size should be a scalar');
    end
    if ~isscalar(order)
        error('The order of AR process should be a scalar');
    end
    if ~isscalar(rho)
        error('The learning rate for step size should be a scalar');
    end
    if ~isscalar(leakage)
        error('The leakage should be a scalar');
    end
    if ~isstring(algo)
        error('The algorithm type should be a scalar');
    end
    if ~isscalar(alpha)
        error('The learning coefficient for Ang & Farhang should be a scalar');
    end
    
    %% Define parameters
    % The order and number of samples
    [~,N] = size(inputSig);
    % The weight
    weight_GASS = zeros(order,N+1);
    % The error
    error_GASS = zeros(1,N);
    % The predicted step size
    mu_GASS = zeros(1,N+1);
    mu_GASS(1) = step0;
    % The phi term
    phi = zeros(order,N+1);
    % The processed signal
    xn_GASS = zeros(order,N);
    
     % Obtain the processed signal x(n-m)
    for m = 1:order
        xn_GASS(m,:) = [zeros(1,m-1),inputSig(1,1:(N-m))];
    end
    
    for n = 1:N
        pred = weight_GASS(:,n).'*xn_GASS(:,n);
        error_GASS(n) = pred - desireSig(n);
        weight_GASS(:,n+1) = (1-mu_GASS(n)*leakage)*weight_GASS(:,n)+mu_GASS(n)*error_GASS(n)*xn_GASS(:,n);
        mu_GASS(n+1) = mu_GASS(n)+rho*error_GASS(n)*xn_GASS(:,n).'*phi(:,n);
        switch algo
            case 'Benveniste'
                phi(:,n+1) = (eye(order)-mu_GASS(n)*xn_GASS(:,n)*xn_GASS(:,n).')*phi(:,n)+error_GASS(n)*xn_GASS(:,n);
            case 'Ang'
                phi(:,n+1) = alpha*phi(:,n)+error_GASS(n)*xn_GASS(:,n);
            case 'Matthews'
                phi(:,n+1) = error_GASS(n)*xn_GASS(:,n);   
        end
    end
    weight_GASS = weight_GASS(:,2:end); 
end