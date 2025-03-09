function [dVdx] = func_grad_V_Jzz(Jz, z)
%FUNC_GRAD_V_JZZ Calculate the gradient of V based on given jacobian and
%transformed states z = phi(x)

% Inputs:
%   Jz              Jacobians of current phi at each datapoint
%                   N-Cellarray, Jz_train{i} is n by n
%
%   z               Transformed states z at which gradient needs to be calculated,
%                   where z(:, i) is the i-th transformed state given by phi(x(:, i))
%                   n by N array


% Outputs:
%   dVdx            Gradients for each datapoint x,
%                   where dVdx(:, i) is the gradient at x(:, i)
%                   n by M array  


    % Get number of states (denoted by N in function description above)
    n_data = size(z, 2);

    % Initialize gradient array
    dVdx = zeros(size(z));

    % Calculate gradient at each datapoint
    for i=1:n_data
        % We always assume the quadratic Lyapunov function 0.1*x.^2
        % Hence, the gradient of V(phi(x)) is given by:
        % dV(dPhi(x))/dx = 0.2 dPhi(x)/dx = 0.2 J_phi * phi(x)
        dVdx(:, i) = 0.2*squeeze(Jz(i, :, :))'*z(:, i);
    end

end

