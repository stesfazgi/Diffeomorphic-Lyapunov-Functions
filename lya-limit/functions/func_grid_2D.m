function [grid_x1, grid_x2, x_grid] = func_grid_2D(x1_lim, x2_lim, n_x1, n_x2)
%FUNC_GRID Helper function to create a 2D grid given per-axis definitions

% Inputs:
%   x1_lim      Array of two elements specifying [x1_min, x1_max] for first
%               axis
%
%   x2_lim      Array of two elements specifying [x2_min, x2_max] for
%               second axis
%
%   n_x1        Integer specifying the number of interpolation points along
%               first axis
%
%   n_x2        Integer specifying the number of interpolation points along
%               second axis


% Outputs:
%   grid_x1     Linspace along first axis
%               1 by n_x1 array
%
%   grid_x2     Linspace along second axis
%               1 by n_x2 array
%
%   x_grid      Grid points
%               2 by n_x1*n_x2 array

    grid_x1 = linspace(x1_lim(1), x1_lim(2), n_x1);
    grid_x2 = linspace(x2_lim(1), x2_lim(2), n_x2);

    [x1_grid, x2_grid] = ndgrid(grid_x1, grid_x2);
    x_grid = [reshape(x1_grid, 1, []); reshape(x2_grid, 1, [])];

end

