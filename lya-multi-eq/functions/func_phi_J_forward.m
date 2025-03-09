function [Jz, z] = func_phi_J_forward(x, phi)
%FUNC_PHI_J_FORWARD Evaluate the jacobian of the full diffeomorphism phi at
%states x

% Inputs:
%   x               States at which to evaluate phi, 
%                   where x(:, i) is the i-th state for evaluation
%                   n by M array
%
%   phi             Current diffeomorphism
%                   See func_phi_init.m for more information


% Outputs:
%   Jz              Jacobians of current phi at each datapoint
%                   M-Cellarray, Jz_train{i} is n by n
%
%   z               Transformed states z = phi(x)
%                   n by M array


    % Get number of layers in diffeo
    n_layers = length(phi.centers_per_layer);

    % Apply complete diffeomorphism for layer 1 to layer "n_layers"
    [Jz, z] = func_phi_J_forward_partial(x, phi, 1, n_layers);

end

