function [alpha_bounds] = func_alpha_bounds_full(H_vec, n_dims, n_centers)
%FUNC_ALPHA_BOUNDS_FULL Calculate the alpha bounds for full covariance matrices

% Inputs:
%   H_vec           Cell array of length "n_centers" consisting 
%                   of "n_dims" by "n_dims" matrices. H_vec{i} is the 
%                   inverse of the covariance for the i-th center
%
%   n_dims          Dimension of state-space
%
%   n_centers       Number of centers


% Outputs:
%   alpha_bounds    N by n matrix containing the elementwise absolute
%                   bounds for each coefficient alpha. For example,
%                   alpha_bounds(i, j) gives the absolute maximum bound
%                   for the coefficent associated with the i-th center
%                   and j-th dimension

    % Initialize bounds
    alpha_bounds = zeros(length(H_vec), n_dims);
    % Iterate over inverse covariances for each center
    for i=1:length(H_vec)
        % Get current inverse covariance
        H = H_vec{i};
        % Perform eigen decomposition
        [Q, D] = eig(H);
        % Get absolute values of eigenvector matrix
        Q_abs = abs(Q);
        % Calculate bounds of diagonalized inverse covariance
        dkdz_bound = exp(-0.5).*sqrt(diag(D));
        % Rotate back
        dkdx_bound = Q_abs*dkdz_bound;
        % Store results in alpha_bounds for current center
        for j=1:n_dims
            alpha_bounds(i, j) = 1/(dkdx_bound(j) * n_centers * n_dims);
        end
    end

end

