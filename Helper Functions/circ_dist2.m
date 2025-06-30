function d = circ_dist2(alpha, beta)
% CIRC_DIST2 Computes circular distance between two sets of angles efficiently
% Usage: d = circ_dist2(alpha, beta)
% Inputs:
%   - alpha: vector of angles (Nx1 or 1xN)
%   - beta:  vector of angles (Mx1 or 1xM)
% Outputs:
%   - d: NxM matrix where d(i,j) = circular difference between alpha(i) and beta(j)

d = angle(exp(1i * (alpha(:) - beta(:)')));  % Efficient circular difference
