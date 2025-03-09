function [loss] = func_fit_lya_mpc_loss_fmincon(Jz_train, z_train, f_target, H_vec, V, alphas)
%FUNC_FIT_LYA_MPC_LOSS_FMINCON Compute the loss for MPC-based diffeomorphic
%Lyapunov search

% Inputs:
%   Jz_train        Jacobians of current phi at each datapoint
%                   N-Cellarray, Jz_train{i} is n by n
%
%   z_train         Datapoints propagated through current phi
%                   n by N array
%
%   f_target        System dynamics at each datapoint (except the origin)
%                   n by N array
%
%   H_vec           Cell array containing the inverse n by n covariance 
%                   matrices for each center/datapoint in x0_train
%                   N-Cellarray, H_vec{i} is n by n
%
%   V               Parameters of Lyapunov function
%                   See func_V_init.m for more information
%
%   alphas          Current guess for the expansion coefficients for all
%                   horizion layers
%                   hor by N by n array


% Outputs:
%   loss            Scalar indicating the loss


    % Get number of horizon steps (denoted by hor in the function description)
    n_hor = size(alphas, 1);
    % Get state space dimension (denoted by n in the function description)
    n_dims = size(alphas, 3);
    % Get number of centers (denoted by N in the function description)
    n_data = size(z_train, 2);

    % Create grammian for kernel evaluation
    kx_data = zeros(n_data, n_data);
    % Create grammian for partial derivatives of kernel evalation
    dkdx_data = zeros(n_data, n_dims, n_data);

    loss = 0;

    % Iterate over horizon steps
    for N=1:n_hor

        % Get alphas for current horizon step
        alpha_curr = squeeze(alphas(N, :, :));
        
        % Calculate grammians
        for j=1:n_data
            H = H_vec{j};
            for i=1:n_data
                d = z_train(:, i) - z_train(:, j);
                s = d'*(H*d);
                kx_data(i, j) = exp(-0.5*s);
                dkdx_data(i, :, j) = -kx_data(i, j)*H*d;
            end
        end

        % Propagate jacobians at each datapoint (=center) through layer of current horizon
        % step, which is parameterized by alphas for current horizon step
        for i=1:n_data
            Jz_train(i, :, :) = (eye(n_dims) + (squeeze(dkdx_data(i, :, :))*alpha_curr)')*squeeze(Jz_train(i, :, :));
        end

        % Propagate datapoints (=centers) through layer of current horizon
        % step, which is parameterized by alphas for current horizon step
        z_train = z_train + (kx_data*alpha_curr)';

    end

    % Calculate normalized gradients
    dVdx_curr = func_grad_V_Jzz(V, Jz_train, z_train);
    dVdx_norm = dVdx_curr ./ vecnorm(dVdx_curr);

    % Calculate Lyapunov constraint violation loss over datapoints
    for i=1:n_data
        loss = loss + dVdx_norm(:, i)'*f_target(:, i);
    end

end

