function [] = polarHeatMap(alpha,r,color)
% Used in Fig. 2
% For polar-planes (alpha,r) plot patches with colors denoted in clr
hold on
for i = 1:size(alpha,1)
    patch_alpha = alpha(i,:); % Pairs of alpha
    patch_alpha = patch_alpha([1,1,2,2,1]); % Coordinates for associated annular
    patch_r = r(i,:); % Pairs of r
    patch_r = patch_r([1,2,2,1,1]); % Coordinates for associated annular
    [x,y] = pol2cart(patch_alpha,patch_r); % Actual annular coordinates
    patch(x,y,color(i,:),'EdgeColor','None');

end
hold off
end