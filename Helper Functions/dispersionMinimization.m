function [min_val,MAD] = dispersionMinimization(input)
% Used in Fig. 3,5 + Extended Fig. 3
% Finds the angle thath minimizes the median absolute deviation (MAD)
% for a set of angular inputs, allowing symmetry (alpha+-pi); not a smooth function, so no point in minimizing by gradient
cands = linspace(0,pi,361*2); % Arbitrary percision
dev = [];
dev = median(min(abs(circ_dist2(input,cands)),abs(circ_dist2(input,cands-pi))),1); % Vectorized
[MAD,idx] = min(dev); % Sometimes a neighborhood around a candidate has equal MAD, choose from among them
if sum(dev == MAD) > 1
    idx = round(median(find((dev-MAD)<10^-10)));
    MAD = dev(idx);
end
min_val = cands(idx);

end