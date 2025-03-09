function [V] = func_Vz(z)
%FUNC_VZ Evaluate the initial guess Lyapunov function V(z) = 0.1*z'*z

% Inputs:
%   z               Transformed states z = phi(x)
%                   n by M array

% Outputs:
%   V               Vector containing the evaluations of V
%                   M by 1 array

    % Get number of centers (denoted by M in the function description)
    n_data = size(z, 2);

    % Evaluate V on all points
    V = zeros(n_data, 1);
    for i=1:n_data
        V(i) = 0.1*sum(z(:, i).^2);
    end

end