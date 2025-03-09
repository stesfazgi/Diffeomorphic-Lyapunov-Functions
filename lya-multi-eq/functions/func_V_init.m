function [V] = func_V_init(v_offset, v_sig, n_centers, v_centers, alphas, H_v)
%FUNC_V_INIT Initialize a M-Lyapunov function constructed from kernels by
%creating a struct with the parameterization

% Inputs:
%   v_offset            Scalar to shift M-Lyapunov function up or down
%
%   v_sig               Covariance scalar for kernels
%
%   n_centers           Number of centers for kernel expansion
%
%   v_centers           Center positions
%                       2 x N array
%
%   alphas              Expansion coefficients
%                       N x 2 array
%
%   H_v                 Inverse of covariance
%                       2 by 2 matrix


% Outputs:
%   V                   Struct containing the parameterization for an
%                       M-Lyapunov function
    

    V.v_offset = v_offset;
    V.v_sig = v_sig;
    V.n_centers = n_centers;
    V.v_centers = v_centers;
    V.alphas = alphas;
    V.H_v = H_v;

end

