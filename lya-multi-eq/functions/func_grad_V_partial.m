function [dVdx] = func_grad_V_partial(V, phi, x, idx_start, idx_end)
%FUNC_GRAD_V_PARTIAL Calculate the partial gradient of V(phi(x)) starting
%and ending at a specific layer. This is useful for example when creating
%morphing videos.

% Examples for phi with 50 layers:
% 
% func_grad_V_partial(V, phi, x, 1, 50) is equal to func_grad_V(phi, x)
%
% func_grad_V_partial(V, phi, x, 1, 30) calculates the gradient of V(phi(x))
% igoring the last 20 layers
%
% func_grad_V_partial(V, phi, x, 20, 50) calculates the gradient of V(phi(x))
% igoring the first 20 layers


% Inputs:
%   V               Parameters of Lyapunov function
%                   See func_V_init.m for more information
%
%   phi             Current diffeomorphism
%                   See func_phi_init.m for more information
%
%   x               States at which to evaluate the gradient, 
%                   where x(:, i) is the i-th state for evaluation
%                   n by M array
%
%   idx_start       Index of layer where to start evaluation of phi
%
%   idx_end         Index of layer where to stop evaluation of phi


% Outputs:
%   dVdx            Gradients for each datapoint x,
%                   where dVdx(:, i) is the gradient at x(:, i)
%                   n by M array

    % Get number of layers in phi
    n_layers = length(phi.centers_per_layer);
    % Get state space dimension (denoted by n in the function description)
    n_dims = size(x, 1);
    % Get number of centers (denoted by M in the function description)
    n_data = size(x, 2);

    % Sanity check
    assert(idx_end <= n_layers, "Requested layer index higher than number of layers");

    % If "up to layer 0" is requested, we perform no transformation
    if idx_end < 1

        % Identity jacobians
        Jx = zeros(n_data, n_dims, n_dims);
        for i=1:n_data
            Jx(i, :, :) = eye(n_dims);
        end

        % Call func_grad_V_Jzz with identity jacobians and z=x (no trafo)
        dVdx = func_grad_V_Jzz(V, Jx, x);

    else

        % Calculate jacobians and transformed states z
        [Jz, z] = func_phi_J_forward_partial(x, phi, idx_start, idx_end);
        % Calculate resulting gradients of V
        dVdx = func_grad_V_Jzz(V, Jz, z);

    end
    


end

