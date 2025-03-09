function [dVdx] = func_grad_V_Jzz(V, Jz, z)
%FUNC_GRAD_V_JZZ Calculate the gradient of V based on given jacobian and
%transformed states z = phi(x)

% Inputs:
%   V               Parameters of Lyapunov function
%                   See func_V_init.m for more information
%
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


    % Get state space dimension (denoted by n in the function description)
    n_dims = size(z, 1);
    % Get number of states (denoted by N in function description above)
    n_data = size(z, 2);

    % Gradient calculation for M-Lyapunov functions is currently only
    % implemented for 2D-State-Spaces
    assert(n_dims == 2, "Only 2D state spaces are supported for this Lyapunov function!")

    % Define the original centers, alphas and inverse covariances of V
    n_centers = V.n_centers;
    v_centers = V.v_centers;
    alphas = V.alphas;
    H_v = V.H_v;

    % Create array for gradients
    dVdx = zeros(n_dims, n_data);

    % Calculate grammian of partial derivative of kernel
    dkdx_data = zeros(n_data, n_dims, n_centers);
    for j=1:n_centers
        for i=1:n_data
            d = (z(:, i) - v_centers(:, j));
            s = d'*(H_v*d);
            dkdx_data(i, :, j) = exp(-0.5*s)*H_v*d;
        end
    end

    % Calculate gradients of V in original space
    for i=1:n_centers
        dVdx = dVdx + squeeze(dkdx_data(:, :, i))'*alphas(i);
    end

    % Project into morphed space by jacobian
    for i=1:n_data
        dVdx(:, i) = squeeze(Jz(i, :, :))'*dVdx(:, i);
    end

end

