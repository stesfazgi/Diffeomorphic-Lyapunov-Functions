function [alpha_star, z_train_H1, loss] = func_fit_lya_mpc_limit(phi, x_train, f_target, H_vec, V, prev_alpha_hor, hor)
%FUNC_FIT_LYA_MPC_LIMIT Runs the optimization to find the expansion coefficients
%"alpha_star" over horizion "hor" to parameterize the diffeomorphism.

% Inputs:
%   phi             Current diffeomorphism
%                   See func_phi_init.m for more information
%
%   x_train         Datapoints (which are also treated as expansion points)
%                   n by N array
%
%   f_target        System dynamics at each datapoint
%                   n by N array
%
%   H_vec           Cell array containing the inverse n by n covariance 
%                   matrices for each center/datapoint in x0_train
%                   N-Cellarray, H_vec{i} is n by n
%
%   V               Parameters of Lyapunov function
%                   See func_V_init.m for more information
%
%   prev_alpha_hor  alpha_star from previous iteration, used for warm start
%                   hor by N by n array
%
%   hor             Integer > 0 indicating the number of horizon steps for
%                   the MPC


% Outputs:
%   alpha_star      Coefficients for each layer in the horizon found by
%                   the solver. alpha_star(1, :, :) are the alphas found
%                   for the first horizon step
%                   hor by N by n array
%
%   z_train_H1      Datapoints x_train propagated through phi, which
%                   are equivalent to the centers used for the layer 
%                   of the first horizon step. These are returned to make
%                   it easier when adding a new layer (which requires
%                   these centers)
%                   n by N array
%
%   loss            Scalar representing the final loss of the optimization

    % Turn off diagnostic messages of fmincon
    options = optimoptions('fmincon','Display','off');
    % Get state space dimension (denoted by n in the function description)
    n_dims = size(x_train, 1);
    % Get number of centers (denoted by N in the function description)
    n_centers = length(H_vec);

    % Create initial alpha guess
    alpha_init = zeros(hor, n_centers, n_dims);

    % Reuse initial guess from previous solution
    for i=2:hor
        alpha_init(i-1, :, :) = prev_alpha_hor(i, :, :);
    end

    % Create alpha bounds.
    % Note: H_vec is kept constant over horizons -> same bounds for all layers
    alpha_lb = zeros(size(alpha_init));
    alpha_ub = zeros(size(alpha_init));
    upper_bound = func_alpha_bounds_full(H_vec, n_dims, n_centers);
    for i=1:hor
        alpha_ub(i, :, :) = upper_bound;
        alpha_lb(i, :, :) = -upper_bound;
    end

    % Calcuate jacobian of forward transformation of current phi
    [Jz_train, z_train_H1] = func_phi_J_forward(x_train, phi);

    % Compile objective expression
    obj = @(a) func_fit_lya_mpc_limit_loss_fmincon(Jz_train, z_train_H1, f_target, H_vec, V, a);

    % Run optimization with box constraints and warm start initial guess
    % Note: all other constraints are included as "soft"-constraints (see func_fit_lya_mpc_limit_loss_fmincon.m)
    [alpha_star, loss] = fmincon(obj, alpha_init, [], [], [], [], alpha_lb, alpha_ub, [], options);

end

