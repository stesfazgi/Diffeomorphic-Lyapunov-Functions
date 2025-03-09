function [phi] = func_phi_init()
%FUNC_PHI_INIT Initialize a diffeomorphism phi by creating a struct with
%empty cell arrays


% Outputs:
%   phi                     Diffeomorphism with 0 layers. The diffeo is
%                           fully specified by the following elements per layer:
%
%   phi.centers_per_layer   Cell-array for the transformed centers of each
%                           layer.
%
%   phi.alphas_per_layer    Cell-array for the coefficients of each layer.
%
%   phi.H_per_layer         Cell-array for the inverse covariances per
%                           layer.
%
%   phi.COV_per_layer       Cell-array for the covariances per layer.


    % Initialize centers per layer. For performance reasons, centers are 
    % always stored in transformed form (they are applied as-is when the 
    % diffeo is evaluated)
    % To get the centers for the i-th layer (assuming state dimension n and
    % number of centers N):
    %   curr_centers = phi.centers_per_layer{i}     (returns n by N array)
    phi.centers_per_layer = {};
    
    % Initialize coefficients per layer.
    % To get the coefficients for the i-th layer (assuming state dimension 
    % n and number of centers N):
    %   curr_alphas = phi.alphas_per_layer{i}       (returns N by n array)
    phi.alphas_per_layer = {};

    % Initialize inverse covariances per layer. For performace reasons, we
    % store the inverse directly to avoid recalculate at every evaluation
    % of phi.
    % To get the inverse covariances for the i-th layer (assuming state
    % dimension n and number of centers N):
    %   curr_H_vec = phi.H_per_layer{i}         (returns Cell-Array of length N of n by n arrays)
    phi.H_per_layer = {};

    % Initialize covariances per layer.
    % To get the covariances for the i-th layer (assuming state
    % dimension n and number of centers N):
    %   curr_H_vec = phi.H_per_layer{i}         (returns Cell-Array of length N of n by n arrays)
    phi.COV_per_layer = {};

end

