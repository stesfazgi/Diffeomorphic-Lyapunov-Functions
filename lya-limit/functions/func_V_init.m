function [V] = func_V_init(v_offset, v_centers, alphas, H_vec)
%FUNC_V_INIT Initialize a M-Lyapunov function constructed from kernels by
%creating a struct with the parameterization

% Inputs:
%   v_offset            Scalar to shift M-Lyapunov function up or down
%
%   v_centers           Center positions
%                       2 x N array
%
%   alphas              Expansion coefficients
%                       N x 2 array
%
%   H_vec               Inverse of covariances per center
%                       N-Cellarray of 2 by 2 matrices


% Outputs:
%   V                   Struct containing the parameterization for an
%                       M-Lyapunov function
    

    V.v_offset = v_offset;
    V.v_centers = v_centers;
    V.alphas = alphas;
    V.H_vec = H_vec;

end

