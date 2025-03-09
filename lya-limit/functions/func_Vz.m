function [V] = func_Vz(V, z)
%FUNC_VZ Evaluate the initial guess Lyapunov function V(z) = 0.1*z'*z

% Inputs:
%   V               Struct containing the parameterization for an
%                   M-Lyapunov function
%
%   z               Transformed states z = phi(x)
%                   n by M array

% Outputs:
%   V               Vector containing the evaluations of V
%                   M by 1 array


    % Get number of datapoints (denoted by M in the function description)
    n_data = size(z, 2);

    % Get number of centers
    n_centers = length(V.H_vec);

    % Evaluate grammian of kernel
    kx_data = zeros(n_data, n_centers);
    for j=1:n_centers
        H_v = V.H_vec{j};
        for i=1:n_data
            d = (z(:, i) - V.v_centers(:, j));
            s = d'*(H_v*d);
            kx_data(i, j) = exp(-0.5*s);
        end
    end

    % Calculate value of V. Note that the kernel is flipped and offset with
    % V.v_offset.
    V = -kx_data*V.alphas + V.v_offset;

end