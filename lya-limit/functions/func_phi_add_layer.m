function [phi] = func_phi_add_layer(phi, z_center, alphas, H_vec, COV_vec)
%FUNC_PHI_ADD_LAYER Add a layer to diffeomorphism phi

% Inputs:
%   phi             Current diffeomorphism
%                   See func_phi_init.m for more information
%
%   z_center        Centers for new layer
%                   Expressed in morhped space: z_center = phi(x_center)
%                   n by N array
%
%   alphas          Expansion coefficients for new layer
%                   N by n array
%
%   H_vec           Cell array containing the inverse n by n covariance 
%                   matrices for each center in z_center
%                   N-Cellarray, H_vec{i} is n by n
%
%   COV_vec         Cell array containing the n by n covariance 
%                   matrices for each center in z_center
%                   N-Cellarray, COV_vec{i} is n by n


% Outputs:
%   phi             New diffeomorphism with added layer


    % Simply extend all cell arrays in the struct with new parameters
    phi.centers_per_layer{end+1} = z_center;
    phi.alphas_per_layer{end+1} = alphas;
    phi.H_per_layer{end+1} = H_vec;
    phi.COV_per_layer{end+1} = COV_vec;

end

