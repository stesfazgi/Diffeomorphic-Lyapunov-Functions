function [Jz, z] = func_phi_J_forward_partial(x, phi, idx_start, idx_end)
%FUNC_PHI_J_FORWARD_PARTIAL Calculate the jacobian of the partial
%diffeomorphism phi between layers idx_start and idx_end at states x


% Inputs:
%   phi             Current diffeomorphism
%                   See func_phi_init.m for more information
%
%   x               States at which to evaluate the jacobian, 
%                   where x(:, i) is the i-th state for evaluation
%                   n by M array
%
%   idx_start       Index of layer where to start evaluation of phi
%
%   idx_end         Index of layer where to stop evaluation of phi


% Outputs:
%   Jz              Jacobians of current phi at each datapoint
%                   M-Cellarray, Jz_train{i} is n by n
%
%   z               Transformed states z = phi(x)
%                   n by M array



    % Get number of layers in diffeo
    n_layers = length(phi.centers_per_layer);
    % Get number of centers (denoted by M in the function description)
    n_data = size(x, 2);
    % Get state space dimension (denoted by n in the function description)
    n_dims = size(x, 1);

    % Sanity check
    assert(idx_end <= n_layers, "Requested layer index higher than number of layers");

    % Initialize jacobians for each state 
    Jz = zeros(n_data, n_dims, n_dims);
    for i=1:n_data
        Jz(i, :, :) = eye(n_dims);
    end

    % Initialization for loop
    z = x;
    
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

        % Initialize grammian for kernel
        kx_data = zeros(n_data, n_centers);
        % Initialize grammian for partial derivative of kernel
        dkdx_data = zeros(n_data, n_dims, n_centers);

        % Calculate grammians for kernel and derivative
        for j=1:n_centers
            H = H_vec{j};
            for i=1:n_data
                d = z(:, i) - z_centers(:, j);
                s = d'*(H*d);
                kx_data(i, j) = exp(-0.5*s);
                dkdx_data(i, :, j) = -kx_data(i, j)*H*d;
            end
        end

        % Update Jacobians at all datapoints.
        % Note: a case destiction is required for n_centers == 1, since
        % matlab does some automatic type conversion...
        if n_centers == 1
            for i=1:n_data
                % Jacobians are propagated by: J_{k+1} = (I + dR/dx)*J_{k}
                % where J_{k} is the jacobian of phi up to layer k
                % dR/dx is the derivative of the residual term
                Jz(i, :, :) = (eye(n_dims) + (dkdx_data(i, :, :)'*alphas)')*squeeze(Jz(i, :, :));
            end
        else
            for i=1:n_data
                Jz(i, :, :) = (eye(n_dims) + (squeeze(dkdx_data(i, :, :))*alphas)')*squeeze(Jz(i, :, :));
            end
        end
        

        % Apply transformation to z by simply adding the residual
        z = z + (kx_data*alphas)';
    end


end

