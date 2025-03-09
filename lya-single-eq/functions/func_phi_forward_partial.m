function [z] = func_phi_forward_partial(z, phi, idx_start, idx_end)
%FUNC_PHI_FORWARD_PARTIAL Evaluate the partial diffeomorphism phi at 
%states x between layers idx_start and idx_end. This is useful for example
%when creating morphing videos


% Inputs:
%   z               States at which to evaluate phi,
%                   where z(:, i) is the i-th state for evaluation
%                   n by M array
%
%   phi             Current diffeomorphism
%                   See func_phi_init.m for more information
%
%   idx_start       Index of layer where to start evaluation of phi
%
%   idx_end         Index of layer where to stop evaluation of phi


% Outputs:
%   z               Transformed states z = phi(x) starting at layer 
%                   idx_start and ending at layer idx_end
%                   n by M array


    % Get number of layers in phi
    n_layers = length(phi.centers_per_layer);
    % Get number of datapoints (denoted by M in the function description)
    n_data = size(z, 2);

    % Sanity check
    assert(idx_end <= n_layers, "Requested layer index higher than number of layers");

    % Iterate over layers in specified index range
    for N=idx_start:idx_end
        % Get centers for current layer
        z_centers = phi.centers_per_layer{N};
        % Get coefficients for current layer
        alphas = phi.alphas_per_layer{N};
        % Get inverse covariances for current layer
        H_vec = phi.H_per_layer{N};
        % Get number of centers for current layer
        n_centers = size(z_centers, 2);
        
        % Calculate grammian
        kx_data = zeros(n_data, n_centers);
        for j=1:n_centers
            H = H_vec{j};
            for i=1:n_data
                d = (z(:, i) - z_centers(:, j));
                s = d'*(H*d);
                kx_data(i, j) = exp(-0.5*s);
            end
        end

        % Apply transformation by simply adding the residual
        z = z + (kx_data*alphas)';

    end

end

