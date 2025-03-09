function [z_center, COV_vec, H_vec] = func_ukf(forward, x_center, x_cov_init, ukf_alpha, ukf_kappa, ukf_beta)
%FUNC_UKF Calcualtes covariance estimates to approx. the morphed kernel by
%using an unscented Kalman filter.

% Inputs:
%   forward             function handle of phi(x)
%
%   x_center            Center points in x space (NOT phi(x))
%                       n by N array
%
%   x_cov_init          Initial covariance matrix for all centers
%                       Note: UKF code currently only supports applying the
%                       same covariance to all centers in x space
%                       n by n array
%
%   ukf_alpha           Scalar, representing the "alpha" param of UKF
%
%   ukf_kappa           Scalar, representing the "kappa" param of UKF
%
%   ukf_beta            Scalar, representing the "beta" param of UKF


% Outputs:
%   z_center            Transformed center points phi(x_center)
%                       n by N array
%
%   COV_vec             Estimated covariances for each center point
%                       N-CellArray of n by n matrices
%
%   H_vec               Inverse of estimated covariances for each center
%                       N-CellArray of n by n matrices


% Additional notes: Transformed mean is phi(z_center), NOT the UKF output.
%                   This guarantees that the centers stay on datapoints,
%                   even if the UKF is parameterized badly
%
%                   This is an implementation of the following paper:
%                   https://groups.seas.harvard.edu/courses/cs281/papers/unscented.pdf



    % Get state space dimension (denoted by n in the function description)
    n_dims = size(x_center, 1);
    % Get number of centers (denoted by N in the function description)
    n_center = size(x_center, 2);

    % Initialize cell arrays to hold results
    COV_vec = {};
    H_vec = {};

    % Perform UKF for each center
    for N=1:n_center

        % Get current center
        x_center_curr = x_center(:, N);
        
        % Calculate "scaling" parameter
        lambda = ukf_alpha^2*(n_dims + ukf_kappa) - n_dims;

        % Calculate UKF sigma points
        x_sig_p = zeros(n_dims, n_dims*2+1);
        x_sig_p(:, 1) = x_center_curr;
        for i=1:n_dims
            XI1 = x_center_curr + (sqrtm((n_dims + lambda)*x_cov_init));
            x_sig_p(:, i+1) = XI1(:, i);
            XI2 = x_center_curr - (sqrtm((n_dims + lambda)*x_cov_init));
            x_sig_p(:, i+1+n_dims) = XI2(:, i);
        end
    
        % Calculate UKF weights
        W0m = lambda/(n_dims+ lambda);
        W0c = W0m + (1 - ukf_alpha^2 + ukf_beta);
        Wi = 1/(2*(n_dims + lambda));

        % Transform sigma points
        z_sig_p = forward(x_sig_p);

        % Calculate UKF estimate of center
        z_ukf_center = W0m*z_sig_p(:, 1);
        for i=1:n_dims*2
            z_ukf_center = z_ukf_center + Wi*z_sig_p(:, i+1);
        end
        
        % Calculate UKF estimate of covariance
        z_cov = W0c*(z_sig_p(:, 1) - z_ukf_center)*(z_sig_p(:, 1) - z_ukf_center)';
        for i=1:n_dims*2
            z_cov = z_cov + Wi*(z_sig_p(:, i+1) - z_ukf_center)*(z_sig_p(:, i+1) - z_ukf_center)';
        end

        % Store result (and also inverse)
        COV_vec{end+1} = z_cov;
        H_vec{end+1} = inv(z_cov);

    end

    % Important: Again, we use phi(x_center) as center, NOT the UKF estimate
    z_center = forward(x_center);

end

