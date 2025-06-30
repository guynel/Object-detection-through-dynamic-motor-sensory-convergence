function [tile_p_norm] = calcPolarProb2D(alpha_r,alpha_edges,r_edges)
% Used in generation of Fig. 2:
% For a 2d (alpha/r) polar visit-rate map, calculate normalized
% probabilities (that take into account annular size)

alpha_vals = alpha_r(:,1); % alpha_values
r_vals = alpha_r(:,2); % r_values


tile_p_norm = [];
for k = 1:(length(r_edges)-1)
    annular_size = pi*((r_edges(k+1).^2)-(r_edges(k).^2));
    idx_r = (r_vals > r_edges(k)) & (r_vals < r_edges(k+1));
    for kk = 1:(length(alpha_edges)-1)
        idx_alpha = (alpha_vals >= alpha_edges(kk)) & (alpha_vals < alpha_edges(kk+1));
        idx = idx_alpha & idx_r;
        tile_p = sum(idx)/length(idx); % Count out of all counts;
        tile_p_norm(k,kk) = tile_p./annular_size; % Normalized to annular size

    end
end
tile_p_norm = tile_p_norm./sum(tile_p_norm(:)); % Scale to sum 1

end