function plot_mad_aligned_ellipse(x, y, scale, N, clr)
% Used in Fig. 1, last panel
% Find for x and y the polar angles, finds the angle for which MAD is
% minimized, project data to both axes of the ellipsed, and plots
% x, y: data vectors
% scale: scale factor for ellipse axes (e.g., 2 for 95% region)
% N: number of ellipse points (default 200)

if nargin < 4
    N = 200;
end
if nargin < 3
    scale = 2;
end

% Stack into 2D points
pts = [x(:), y(:)];

% Get angles from origin
[alpha, ~] = cart2pol(pts(:,1), pts(:,2));

% Find the direction with minimal angular MAD
theta = dispersionMinimization(alpha);

% Compute angular deviations from the two axes
d_primary = median(min(abs(circ_dist2(alpha,theta)),abs(circ_dist2(alpha,theta-pi))),1); % Vectorized
theta_orth = theta+pi/2;
d_ortho = median(min(abs(circ_dist2(alpha,theta_orth)),abs(circ_dist2(alpha,theta_orth-pi))),1); % Vectorized

% Compute MAD of angular deviations
a = scale * d_primary;
b = scale * d_ortho;

% Parametrize ellipse in aligned frame
t = linspace(0, 2*pi, N);
% Construct ellipse aligned to major axis (horizontal)
ellipse = [a*cos(t); b*sin(t)];

% Now rotate the ellipse *into the original frame*, by +theta
R = [cos(theta-pi/2), -sin(theta-pi/2); 
     sin(theta-pi/2),  cos(theta-pi/2)];
ellipse_xy = R * ellipse;

% Plot (optional)
hold on
patch(ellipse_xy(1,:), ellipse_xy(2,:), clr, 'FaceAlpha', 0.4);
plot(ellipse_xy(1,:), ellipse_xy(2,:), 'Color', clr, 'LineWidth', 2);
hold off;
end
