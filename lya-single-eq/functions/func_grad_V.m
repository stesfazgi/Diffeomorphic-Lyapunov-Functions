function [dVdx] = func_grad_V(phi, x)
%FUNC_GRAD_V Calculate the gradient of the diffeomorphic Lyapunov function
%at states x. ( dVdx = d V(phi(x)) / d x )

% Inputs:
%   phi             Current diffeomorphism
%                   See func_phi_init.m for more information
%
%   x               States at which to evaluate the gradient, 
%                   where x(:, i) is the i-th state for evaluation
%                   n by M array


% Outputs:
%   dVdx            Gradients for each datapoint x,
%                   where dVdx(:, i) is the gradient at x(:, i)
%                   n by M array


    % Get number of layers in diffeo
    n_layers = length(phi.centers_per_layer);
    
    % Call partial gradient calculation with full layer count (1 to n_layers)
    % to get the gradient of the complete diffeo
    dVdx = func_grad_V_partial(phi, x, 1, n_layers);

end

