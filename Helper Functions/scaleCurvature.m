function [allJobs,scale_coeff] = scaleCurvature(allJobs)
% Used after Fig. 2
% Scales curvature such that the overall distribution of alpha, derived
% from theta-dot and kappa-dot is as uniform as possible
stats = [];
% Gather all paired observations
for trial_num = 1:length(allJobs)
    d_ang = allJobs{trial_num}.Tracks.diffAngles(:);
    d_curv = allJobs{trial_num}.Tracks.diffCurv(:);
    nan_idx = isnan(d_ang) | isnan(d_curv);

    d_ang(nan_idx) = [];
    d_curv(nan_idx) = [];

    stats = [stats; [d_ang,d_curv]];
end

% Find scaling factor
res = [];
scale_vals = 0.01:0.01:2; % Hand-picked, could also be minimized, doesn't change reuslt
for c = scale_vals
    theta = atan2(c*stats(:,2), stats(:,1));
    empirical_cdf = histcounts(theta,-pi:0.1:pi,'Normalization','cdf');
    unif_cdf = linspace(0,1,length(empirical_cdf));
    x_steps = linspace(-pi,pi,length(empirical_cdf));
    res(end+1) = trapz(x_steps,abs(unif_cdf-empirical_cdf));
end
[~,idx] = min(res);
scale_coeff = scale_vals(idx); % Scale factor

% Apply to kappa-dot everywhere
for trial_num = 1:length(allJobs)
    d_curv = allJobs{trial_num}.Tracks.diffCurv;
    d_curv = scale_coeff.*d_curv;
    allJobs{trial_num}.Tracks.diffCurv = d_curv;

end

end