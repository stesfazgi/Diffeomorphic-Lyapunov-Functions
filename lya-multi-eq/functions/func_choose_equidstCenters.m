function [x_centers] = func_choose_equidstCenters(x_data, n_centers)
%FUNC_CHOOSE_EQUIDSTCENTERS Helper function to choose approx. equally
%spaced datapoints from trajectory

% Inputs:
%   x_data          Trajectory of N datapoints in n dimensional statespace
%                   n by N array
%
%   n_centers       Number of points to pick from trajectory


% Outputs:
%   x_centers       Approx. equally spaced datapoints chosen from trajectory
%                   n by n_centers array

    x_data = unique(x_data.','rows','stable').';

    n_dims = size(x_data, 1);
    n_data = size(x_data, 2);

    t = [0];
    x_data_unique = [x_data(:, 1)];
    
    for i=2:n_data
        t_new = norm(x_data(:, i-1) - x_data(:, i), 2);
        if t_new > 1e-3
            t_new = t_new + t(end);
            t = [t t_new];
            x_data_unique = [x_data_unique x_data(:, i)];
        end
        
    end

    if size(x_data_unique, 2) < n_centers
        warning("Not enough unique states to choose centers! Filling with last state");
        n_unique = size(x_data_unique, 2);
        x_centers = repmat(x_data(:, end), 1, n_centers);
        x_centers(:, 1:n_unique) = x_data_unique;
        
    else
        x_centers = zeros(n_dims, n_centers);
        for i=1:n_dims
            x_centers(i, :) = interp1(t, x_data_unique(i, :), linspace(t(end), 0, n_centers));
        end
    end

end

