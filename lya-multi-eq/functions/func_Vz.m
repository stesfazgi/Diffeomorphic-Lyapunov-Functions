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


    % Get number of centers (denoted by M in the function description)
    n_data = size(z, 2);

    % Evaluate grammian of kernel
    kx_data = zeros(n_data, V.n_centers);
    for j=1:V.n_centers
        for i=1:n_data
            d = (z(:, i) - V.v_centers(:, j));
            s = d'*(V.H_v*d);
            kx_data(i, j) = exp(-0.5*s);
        end
    end

    % Calculate value of V. Note that the kernel is flipped and offset with
    % V.v_offset.
    V = -kx_data*V.alphas + V.v_offset;

end