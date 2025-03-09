function [] = func_run_seq(params, output_dir)
%FUNC_RUN_SEQ Searches for a diffeo Lyapunov function, saves the result
%to a .mat file and creates a .jpg figure with the results. This is the
%script version of the notebook "run_single_diffeo_lya.mlx" to make it
%easier to run multiple experiments in a row

% Inputs:
%   params.shape_id         ID of LASA shape
%   params.n_hor            Number of MPC horizon steps
%   params.sig              Sigma scaling factor for diffeo covariances
%
%   output_dir              String specifying the output directory for
%                           the .mat and .jpg files


% Outputs:
%   -                       (The script saves all results to files)

    % Get paramters for current run
    shape_id = params.shape_id;
    n_hor = params.n_hor;
    sig = params.sig;

    % Load dataset
    n_samples = 30;
    n_data = n_samples;
    [data_pos, data_vel, shapename, ~] = plot_shape(shape_id, 1, n_samples, false, [50 50]);
    % Normalize dataset
    scale_factor = 1/max(abs(data_pos(:)));
    data_pos = data_pos*scale_factor;
    data_vel = data_vel*scale_factor;
    x_data = data_pos;
    dx_data = data_vel;

    % Params for phi
    n_iterations = 50;
    ukf_alpha = 1;
    ukf_kappa = 1;
    ukf_beta = 2;
    cov_init = eye(2)*sig;
    H_init = diag(1./diag(cov_init));
    
    % Initialize diffeo
    phi = func_phi_init();
    
    % Initialize alpha horizon
    alpha_star = zeros(n_hor, n_data, 2);

    % Initialize logging arrays
    loss_per_iter = zeros(n_iterations, 1);
    vios_per_iter = zeros(n_iterations, 1);
    
    % Normalize system velocities
    dx_norm_data = dx_data ./ vecnorm(dx_data);
    
    tic;
    for N=1:n_iterations
    
        % Get transformed covariances
        forward = @(x) func_phi_forward(x, phi);
        [~, COV_vec, H_vec] = func_ukf(forward, x_data, cov_init, ukf_alpha, ukf_kappa, ukf_beta);
    
        % Run MPC
        [alpha_star, z_train_H1, loss] = func_fit_lya_mpc(phi, x_data, dx_norm_data, H_vec, alpha_star, n_hor);
    
        % Add new layer (first horizon step of MPC)
        phi = func_phi_add_layer(phi, z_train_H1, squeeze(alpha_star(1, :, :)), H_vec, COV_vec);
    
        % Check number of violations
        n_vio = 0;
        dVdx_curr = func_grad_V(phi, x_data);
        for i=1:n_data-1
            L2_curr = dVdx_curr(:, i)'*dx_data(:, i);
            if L2_curr > 0
                n_vio = n_vio + 1;
            end
        end
    
        % Logging
        vios_per_iter(N) = n_vio;
        loss_per_iter(N) = loss;
        fprintf("ITER=%i, Vios=%i, Loss=%f \n", N, n_vio, loss);
    
    end
    t_opt = toc;

    % Bounds for plotting
    x1_lim = [-1.5; 1.5];
    x2_lim = [-1.5; 1.5];
    
    % Number of evaluation points per dimension
    n_eval_x1 = 100;
    n_eval_x2 = 100;
    
    % Create evaluation grid
    [evalgrid_x1, evalgrid_x2, x_eval] = func_grid_2D(x1_lim, x2_lim, n_eval_x1, n_eval_x2);
    n_eval = size(x_eval, 2);

    % Evaluate diffeo on eval grid
    z_eval = func_phi_forward(x_eval, phi);
    
    % Evaluate V on transformed evaluation points
    V_fit = func_Vz(z_eval);
    V_fit = reshape(V_fit, size(evalgrid_x1, 2), size(evalgrid_x2, 2));
    
    % Plot the result
    fig = figure("Position", [0 0 3*700 700]);
    subplot(1, 3, 1);
    hold on;
    plot(x_data(1, :), x_data(2, :), "rx", "LineStyle", "none");
    contour(evalgrid_x1, evalgrid_x2, V_fit', linspace(0, max(V_fit(:)), 50));
    colorbar;
    subplot(1, 3, 2);
    hold on;
    plot(linspace(1, n_iterations, n_iterations), loss_per_iter);
    xlabel("Iteration");
    ylabel("Loss");
    subplot(1, 3, 3);
    hold on;
    plot(linspace(1, n_iterations, n_iterations), vios_per_iter);
    xlabel("Iteration");
    ylabel("Violations");
    sgtitle(sprintf("Shape=%s, NVios=%i, Loss=%.1f, sig=%.2f", shapename, vios_per_iter(n_iterations), loss_per_iter(n_iterations), sig));

    % Create paths for .mat and .jpg files
    fig_path = sprintf("%s/%s.jpg", output_dir, shapename);
    mat_path = sprintf("%s/%s.mat", output_dir, shapename);

    % Save plot to .jpg
    saveas(fig, fig_path);
    % Save result to .mat
    save(mat_path, "phi", "x_data", "dx_data", "params", "t_opt", "vios_per_iter", "loss_per_iter", "n_hor", "n_iterations", "sig");
    close all;

end

