function [allJobs] = defineApproachPhase(allJobs,px2mm,cont_thresh_px,extend_to_cycle_end_flag)
%% Sandbox to test cutoffs for approach-phase
snout_obj_cont_thresh_px = cont_thresh_px; % Threshold that defines snout-object contact. Divide this by px2mm to get quantity in mm

d_step = 1;
head_center_y = cellfun(@(x) x.Tracks.headCenter(2,:),allJobs,'Uni',0);
warning('off')
lp_head_center_y = cellfun(@(x) lowpassFilt(x,'Fst',100),head_center_y,'Uni',0);
d_head_center_y = cellfun(@(x) multiStepDiff(x,d_step),lp_head_center_y,'Uni',0);
warning('on')

% Find minimal head-velocity from head-center, data-driven

[y,x] = histcounts(log(abs(cell2mat(d_head_center_y))),'Normalization','pdf');
opts = statset('MaxIter', 1000);
warning('off')
gm = fitgmdist(log(abs(cell2mat(d_head_center_y)))',3,'Options',opts);
warning('on')
x_fine = linspace(-8,3, 500);
% Evaluate the GMM
pdf_gm = pdf(gm, x_fine(:));
% Scale the GMM PDF to match the histogram area
scaled_pdf_gm = pdf_gm * sum(diff(x)) * sum(y);
% Visualize fit:
%{
figure
histogram(log(abs(cell2mat(d_head_center_y))),'Normalization','pdf')
yyaxis right
hold on
% Plot the GMM on top of the histogram
plot(x_fine, scaled_pdf_gm, 'b-', 'LineWidth', 2);
% Plot individual components of the GMM
for k = 1:gm.NumComponents
component_pdf = gm.ComponentProportion(k) * normpdf(x_fine, gm.mu(k), sqrt(gm.Sigma(k)));
scaled_component_pdf = component_pdf * sum(diff(x)) * sum(y);
plot(x_fine, scaled_component_pdf, 'r-','LineWidth',2);
end
set(gca,'YTick',[])
yyaxis left;
xlabel('Y-axis Head-Velocity ');
ylabel('Frequency');
title('Multi-Mixture Guassian Fit');
legend('Histogram', 'GMM', 'Components');
%}
% Assign as the mode of the left-most Gaussian
stopping_thresh_px = exp(min(gm.mu));
min_head_y_vel_px = stopping_thresh_px;

% Find head-pauses
head_pause = cellfun(@(x) find(abs(x)<min_head_y_vel_px,1),d_head_center_y,'Uni',0);
% Find snout-object contacts
first_snout_cont = cellfun(@(x) find(x.Tracks.snoutObjDist < ...
    snout_obj_cont_thresh_px,1),allJobs,'Uni',0);


for k = 1:length(head_pause)
    if isempty(head_pause{k})
        head_pause{k} = NaN;
    end
    if isempty(first_snout_cont{k})
        first_snout_cont{k} = NaN;
    end
end
head_pause = cell2mat(head_pause);
first_snout_cont = cell2mat(first_snout_cont);

% Assigne end of approach (called trim_idx)
trim_idx = first_snout_cont; % Assign as default snout-object contact, if occurs
trim_idx(isnan(trim_idx)) = head_pause(isnan(trim_idx)); % Otherwise assign head-pause

% Initially assign this definition
for k = 1:length(allJobs)
    allJobs{k}.Tracks.trim_idx = trim_idx(k);
    allJobs{k}.Tracks.trim_by_snout_cont = ~isnan(first_snout_cont(k));
end


% Get cycle division (depends on trim
if ~isfield(allJobs{1}.Tracks,'cycles')
     [~,allJobs] = extractWhiskerCycles(allJobs,cont_thresh_px);
end


end