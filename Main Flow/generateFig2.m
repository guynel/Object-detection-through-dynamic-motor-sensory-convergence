function [allJobs,cont_thresh,alpha_star_val] = generateFig2(allJobs,px2mm,curv_scale_coeff)
get_permutation_p_val_flag = 1; % Do you want to do bootstrap test? Resoure-heavy...
%% Parameters
% First visit-rate map is pre-cont or all?
pre_cont_flag = true;

% log or linear colors
log_scale_clrs = true;

% do we include nose contact?
excludeNoseTipCont = true;
noseTipContThresh = 0.1; % What's the minimal distance sufficient to be called a nose contact, in cm

% Pre-contact threshold
%pre_cont_thresh = 0.5;
pre_cont_thresh = 0.134; % This actually is calculated empirically below, it's the threshold in Fig. 2a, in units of cm; the output cont_thresh is this value 

%% Get some data
rawWhObjDist = cell2mat(cellfun(@(x) x.Tracks.whiskerObjDist(:)',allJobs,'Uni',0));
rawWhObjDist = rawWhObjDist./(px2mm*10);
snObjDist = cell2mat(cellfun(@(x) x.Tracks.snoutObjDist(:)',allJobs,'Uni',0));
snObjDist = cellfun(@(x) x.Tracks.snoutObjDist,allJobs,'Uni',0);
n_row = cellfun(@(x) size(x.Tracks.whiskerObjDist,1),allJobs,'Uni',0);
snObjDist = cellfun(@(x,y) repmat(x,y,1),snObjDist,n_row,'Uni',0);
snObjDist = cell2mat(cellfun(@(x) x(:)',snObjDist,'Uni',0));
snObjDist = snObjDist./(px2mm*10);
% Remove frame with snout-object contact
if excludeNoseTipCont
    rawWhObjDist(snObjDist<noseTipContThresh) = NaN;
end

pre_cont_dist = cellfun(@(x) repmat((cumsum(any((x.Tracks.whiskerObjDist./(px2mm*10)) < pre_cont_thresh, 1), 2) == 0),size(x.Tracks.whiskerObjDist,1),1), allJobs, 'Uni',0);
pre_cont_dist = cell2mat(cellfun(@(x) x(:)',pre_cont_dist,'Uni',0));

% Get derivatives and distance, and remove NaNs
d_curv = cell2mat(cellfun(@(x) x.Tracks.diffCurv(:)',allJobs,'Uni',0));
d_ang = cell2mat(cellfun(@(x) x.Tracks.diffAngles(:)',allJobs,'Uni',0));
nan_vals = isnan(rawWhObjDist) | isnan(d_ang) | isnan(d_curv);
d_curv = d_curv(~nan_vals);
d_ang = d_ang(~nan_vals);
whObjDist = rawWhObjDist(~nan_vals);
pre_cont_dist = pre_cont_dist(~nan_vals);


%% Make figure
figure_panel = figure;
pos = [0 0 21 29.7];
set(figure_panel, 'PaperUnits','centimeters');
set(figure_panel,'Units','Centimeters','Position',pos);
set(figure_panel, 'PaperSize', [pos(3) pos(4)]);
set(figure_panel,'PaperPosition',pos)
figure_handles = gobjects(0);
figure_handles(1) = figure_panel;

%% Distribution of distance, which will determine contact
fh(1) = axes();
hold on
plot(NaN,NaN,'r','LineWidth',1.5)
plot(NaN,NaN,'b','LineWidth',1.5)
xx = -10:0.01:2.5;
xx_bar = -10:0.25:2.5;
yy = histcounts(log(whObjDist),xx_bar);
yy = yy./(sum(yy)*25);
b_bar = bar(xx_bar(1:end-1)+0.5*mode(diff(xx_bar)),yy);
set(b_bar,'FaceColor',[0 0 0],'FaceAlpha',0.1,'EdgeColor','None')
hold on
yy = ksdensity(log(whObjDist),xx,'Bandwidth',0.15);
yy = yy./sum(yy);
plot(xx,yy,'Color',[0 0 0 0.2],'LineWidth',2)
dyy = multiStepDiff(yy,1);
distrib_idx = find(dyy>0 & [NaN,dyy(1:end-1)]<0,1,'last');
plot(xx(1:distrib_idx),yy(1:distrib_idx),'Color','r','LineWidth',2)
plot(xx(distrib_idx:end),yy(distrib_idx:end),'Color','b','LineWidth',2)
plot(xx(distrib_idx)*[1 1],max(ylim)*[0,1],'k--')
fh(1).Units = 'Centimeters';
fh(1).Position = [2.25 23.5 3 1];
xlabel('Dist._{(Whisker, Object)} [cm]');
ylabel('Prob.');
xlim([-4 2.5])
ylim([10^-4,10^-2])
set(gca,'YScale','log');
cont_thresh = exp(xx(distrib_idx));

% Everything is done in natural log, but we want units ticks to be in log10
% of cm, so this is a slight work around for that (base of transformation
% doesn't make a difference either way):
x_ticks = 10.^[-1:2];
x_tick_lbls = arrayfun(@(x) x, x_ticks, 'Uni',0);
set(gca,'XTick',log(x_ticks),'XTickLabel',x_tick_lbls)

leg = legend({'Contact';'Non-contact'});
set(leg,'Units','centimeters')
set(leg,'Position',[2.4 24.75 2.725 0.8])

%% Parameters for 2D plots

n_polar_wedges = 65; % Theta ticks 
n_r_circles = 31; % r ticks
num_prob_categories_joint = 20; % Colormap bins
max_r_val_diff = 3; % maximal r - also determines marginals
alpha_joint_plot_range = linspace(-pi,pi,n_polar_wedges); % alpha bins
r_joint_plot_range = linspace(0,max_r_val_diff,n_r_circles); % radius bins
ang_marginal_plot_range = -max_r_val_diff:0.15:max_r_val_diff; % d_ang marginals
curv_marginal_plot_range = -max_r_val_diff:0.15:max_r_val_diff; % d_curv marginals
prob_clrs = gray(num_prob_categories_joint); % Colormap
if log_scale_clrs
    prob_bins = linspace(-4,-2.75,num_prob_categories_joint);
    clr_ticks = -4:0.5:-2.75;
else
    prob_bins = linspace(0,0.0015,num_prob_categories_joint);
end


%% 2D visit-rates: All, cont and non-cont
figure(figure_panel);
main_fig_overall_Ns = [];

% Gather count data
counts_for_heatmap = cell(1,3);
for mm = 1:3 % Iterate: pre-cont, cont, no-cont
    [alpha_vals,r_vals] = cart2pol(d_ang,d_curv);
    
if mm == 1
    if pre_cont_flag
        alpha_vals(~pre_cont_dist) = [];
        r_vals(~pre_cont_dist) = [];
    end
elseif mm == 2
    alpha_vals(pre_cont_dist | (whObjDist>cont_thresh)) = [];
    r_vals(pre_cont_dist | (whObjDist>cont_thresh)) = [];
elseif mm == 3
    alpha_vals(pre_cont_dist | (whObjDist<=cont_thresh)) = [];
    r_vals(pre_cont_dist | (whObjDist<=cont_thresh)) = [];
end
main_fig_overall_Ns(end+1) = length(alpha_vals);
alpha_vals_binned = discretize(alpha_vals,alpha_joint_plot_range);
r_val_binned = discretize(r_vals,r_joint_plot_range);

counts_for_heatmap{mm} = NaN(n_polar_wedges-2,n_r_circles-2);
for k = 1:(length(alpha_joint_plot_range)-1)
    for kk = 1:(length(r_joint_plot_range)-1)
        counts_for_heatmap{mm}(k,kk) = sum(alpha_vals_binned == k & r_val_binned == kk);
    end
end
end

% Get probabilities from counts:
probs_for_heatmap = cellfun(@(x) x./nansum(x(:)),counts_for_heatmap,'Uni',0);

for mm = 1:3 % Iterate again for plotting
fh(mm+1) = axes();
for k = 1:(length(alpha_joint_plot_range)-1)
    for kk = 1:(length(r_joint_plot_range)-1)

        x = alpha_joint_plot_range([k,k+1,k+1,k,k]);
        y = r_joint_plot_range([kk,kk,kk+1,kk+1,kk]);
        [x,y] = pol2cart(x,y);
        if isnan(probs_for_heatmap{mm}(k,kk))
            clr_idx = 1;
        else
            if log_scale_clrs
                [~,clr_idx] = min(abs(log10(probs_for_heatmap{mm}(k,kk))-prob_bins));
            else
                [~,clr_idx] = min(abs(probs_for_heatmap{mm}(k,kk)-prob_bins));
            end
        end
        patch(x,y,prob_clrs(clr_idx,:),'EdgeColor','None');
    end
end

        
% Add colored borders
border_circ_alpha = linspace(-pi,pi,37);
border_circ_alpha = [border_circ_alpha(1:end-1);border_circ_alpha(2:end)]';
border_circ_r = (max_r_val_diff.*[1 1]+[0 0.2]).*ones(size(border_circ_alpha));
if mm == 1
    circ_clr = [0 0 1];
elseif mm == 2
    circ_clr = [1 0 0];
elseif mm == 3
    circ_clr = [0 0 1];
end
border_clrs = repmat(circ_clr,size(border_circ_alpha,1),1);
polarHeatMap(border_circ_alpha,border_circ_r,border_clrs);


% Add alpha-star line
if mm == 2
    pos_alpha = alpha_vals;
    pos_alpha(pos_alpha<0) = pos_alpha(pos_alpha<0)+pi;
    [y,x] = ksdensity(pos_alpha,0:0.001:pi);
    [~,max_idx] = max(y);
    max_val = x(max_idx);
    yy = -max_r_val_diff:0.01:max_r_val_diff;
    hold on
    plot(yy./tan(max_val),yy,'r--','LineWidth',1.5);
    alpha_star_val = max_val;
end

% Adjust axes, ticks, labels
fh(mm+1).Units = 'Centimeters';
if mm ~=1
    set(fh(mm+1),'TickLength',[0 0],'xticklabels',{},'yticklabels',{});
    set(gca,'Color','None','XColor','None','YColor','None');
else
    set(gca,'Color','None','XColor','k','YColor','k');
end
axis([-1 1 -1 1]*(max_r_val_diff+0.2))
axis equal
if mm==1

    colormap(prob_clrs);
    clr_bar = colorbar();
    clr_bin_centers = linspace(0,1,1+length(prob_clrs));
    clr_bin_centers = clr_bin_centers(1:end-1)+0.5*mode(diff(clr_bin_centers));
    if log_scale_clrs
        y_ticks = interp1(10.^prob_bins, clr_bin_centers, 10.^(clr_ticks), 'linear', 'extrap');
        y_lbls =  arrayfun(@(x) ['10^{',num2str(x),'}'],clr_ticks,'Uni',0);
    else
        y_ticks = [0,1];
        y_lbls = arrayfun(@(x) num2str(x),[min(prob_bins),max(prob_bins)],'Uni',0);
    end
    set(clr_bar,'YTick',y_ticks,'YTickLabel',y_lbls);
    set(clr_bar,'Units','centimeters','Position',[5.5 17.65 0.2 3.5]);
    clr_bar.Label.String = 'Probability Mass';

end
end

% Adjust positions
fh(2).InnerPosition = [1.75 17.625 3.5 3.5];
axes(fh(2))
axis([-1 1 -1 1]*(max_r_val_diff+0.2))
set(gca,'XColor','k','YColor','k');
x_ticks = xticks;
x_tick_labels = x_ticks*500/1000;
set(gca,'XTick',x_ticks,'XTickLabel',x_tick_labels)
y_labels_display = (yticks/curv_scale_coeff)*500/1000;
scale_log_ten = floor(log10(max(abs(y_labels_display))));
y_labels_rounded = round(y_labels_display ./ (10^scale_log_ten)) * (10^scale_log_ten);
y_ticks_new = (y_labels_rounded * 1000 / 500) * curv_scale_coeff;
set(gca, 'YTick', y_ticks_new, 'YTickLabel', y_labels_rounded)
ylabel('$\dot{\kappa}$ $[m^{-1}/s]$ $\times$ $10^3$','rotation',90,'Interpreter','latex')
xlabel('$\dot{\theta}$ $[degrees/s]$  $\times$ $10^3$','Interpreter','latex')

fh(3).InnerPosition = [8 20.5 3.5 3.5];
axes(fh(3))
axis([-1 1 -1 1]*(max_r_val_diff+0.2))
fh(4).InnerPosition = [8 14.75 3.5 3.5];
axes(fh(4))
axis([-1 1 -1 1]*(max_r_val_diff+0.2))

txt1 = text(fh(2), 0.5, 1.05, 'Pre-contact', 'Units', 'normalized', ...
     'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', 'FontSize', 12);
txt2 = text(fh(3), 0.5, 1.6, 'Contact', 'Units', 'normalized', ...
     'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', 'FontSize', 12);
txt3 = text(fh(4), 0.5, 1.05, {'Non-contact';'(after 1st cont.)'}, 'Units', 'normalized', ...
     'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', 'FontSize', 12);

%% Statistics test for the differences in the 2D distributions
rng(121987,'twister') % Set random seed for reproducibility

% Gather data in one big matrix
cont_labels_xy = cell(1,3);
for mm = 1:3
    [alpha_vals,r_vals] = cart2pol(d_ang,d_curv);

    if mm == 1
        if pre_cont_flag
            alpha_vals(~pre_cont_dist) = [];
            r_vals(~pre_cont_dist) = [];
        end
        
    elseif mm == 2
        alpha_vals(pre_cont_dist | (whObjDist>cont_thresh)) = [];
        r_vals(pre_cont_dist | (whObjDist>cont_thresh)) = [];
    elseif mm == 3
        alpha_vals(pre_cont_dist | (whObjDist<=cont_thresh)) = [];
        r_vals(pre_cont_dist | (whObjDist<=cont_thresh)) = [];
    end
    cont_labels_xy{mm} = [alpha_vals(:),r_vals(:)];
end

% Add labels for pre-cont, cont, no-cont
grpd_cont_labels = [];
for mm = 1:length(cont_labels_xy)
    grpd_cont_labels = [grpd_cont_labels;[cont_labels_xy{mm},mm*ones(size(cont_labels_xy{mm},1),1)]];
end

% Generate probability distributions in 2d
probMaps = cell(1,3);
for k = 1:3
    idx = grpd_cont_labels(:,3) == k;
    probMaps{k} = calcPolarProb2D(grpd_cont_labels(idx,1:2),alpha_joint_plot_range,r_joint_plot_range);
end
% Function for Jensen-Shannon Divergence:
JSdiv = @(P,Q) 0.5 * (nansum(P .* log2(P ./ ((P + Q)/2 ))) + nansum(Q .* log2(Q ./ ((P + Q)/2 ))));
% JSD(cont,no-cont)-JSD(no-cont,pre-cont)
gt = JSdiv(probMaps{2}(:),probMaps{3}(:))-JSdiv(probMaps{1}(:),probMaps{3}(:));

% Perform a permutation test, is our ground-truth (gt) more than expected
% for shuffled data?
sim_res = [];
if get_permutation_p_val_flag
for iter = 1:1000
    probMaps = cell(1,3);
    shuff_labels = grpd_cont_labels;
    shuff_labels(:,3) = shuff_labels(randperm(length(shuff_labels(:,3))),3);
    for k = 1:3
        idx = shuff_labels(:,3) == k;
        probMaps{k} = calcPolarProb2D(shuff_labels(idx,1:2),alpha_joint_plot_range,r_joint_plot_range);
    end
    sim_res(iter) = JSdiv(probMaps{2}(:),probMaps{3}(:))-JSdiv(probMaps{1}(:),probMaps{3}(:));
end
p_val_JSD_diffs = mean(sim_res>gt);

end



%% Marginal for theta-dot and kappa-dot
% Here we force the conversion to units\s by multiplying everything by 500,
% but it doesn't really matter in any of the ratio-based calculations
line_clr =  'rb';
marg_fh = gobjects(0);
% Side marginal, kappa distributions for cont/no-cont
for k = 1:2
    marg_fh(end+1) = axes();
    if k == 1
        [y,x] = histcounts(500*d_curv(whObjDist<cont_thresh),500*curv_marginal_plot_range,'Normalization','pdf');
    else
        [y,x] = histcounts(500*d_curv(whObjDist>cont_thresh),500*curv_marginal_plot_range,'Normalization','pdf');
    end
    y = y./sum(y(:));
    x = x(1:end-1)+0.5*mode(diff(x));
    plot(y,x,'Color',line_clr(k));
    hold on
    ylim([min(500*curv_marginal_plot_range),max(500*curv_marginal_plot_range)])
    set(marg_fh(end),'yaxislocation','right');
    joint_y_lim = get(fh(3+(k-1)),'YLim');
    joint_cm =  get(fh(3+(k-1)),'Position');
    scale = (2 * joint_y_lim(2)) / joint_cm(end);                  % data units per cm
    new_height_cm = (2 *  max(r_joint_plot_range)) / scale;           % height needed for same visual scale
    center_y = joint_cm(2) + joint_cm(end) / 2;                         % vertical center of first axis
    new_y0 = center_y - new_height_cm / 2;                  % align centers
    if k == 2
        set(marg_fh(end),'xaxislocation','top','XTickLabelRotation',0);
    end
    marg_fh(end).Units = 'Centimeters';
    marg_fh(end).Position = [11.5, new_y0, 1, new_height_cm];
    xlabel('Prob.');
    y_lbl = ylabel('$\dot{\kappa}$ $[m^{-1}/s]$ $\times$ $10^3$','rotation',90,'Interpreter','latex');
    set(y_lbl,'Units','centimeters');
    set(y_lbl,'Position',get(y_lbl,'Position')+[-0.4 0 0 ])
    set(marg_fh(end),'YTickLabelRotation',-45);

    new_y_ticks = round(get(gca,'YTick')./(curv_scale_coeff*1000))*1000;
    set(gca,'YTick',new_y_ticks.*curv_scale_coeff);
    set(gca,'YTickLabel',new_y_ticks/1000);
end
max_x_lim = max([marg_fh(end-1:end).XLim]);
xlim(marg_fh(end-1:end),[0,max_x_lim]);

% Top marginal, angle distributions for cont/no-cont
line_clr =  'rb';
for k = 1:2
    marg_fh(end+1) = axes();
    if k == 1
        [y,x] = histcounts(500*d_ang(whObjDist<=cont_thresh),500*ang_marginal_plot_range,'Normalization','pdf');
    else
        [y,x] = histcounts(500*d_ang(whObjDist>cont_thresh),500*ang_marginal_plot_range,'Normalization','pdf');
    end
    y = y./sum(y(:));
    x = x(1:end-1)+0.5*mode(diff(x));
    plot(x,y,'Color',line_clr(k));
    hold on
    xlim([min(500*ang_marginal_plot_range),max(500*ang_marginal_plot_range)])
    joint_x_lim = get(fh(3+(k-1)),'XLim');
    joint_cm =  get(fh(3+(k-1)),'Position');
    scale = (2 * joint_x_lim(2)) / joint_cm(3);                  % data units per cm
    new_width_cm = (2 *  max(r_joint_plot_range)) / scale;           % height needed for same visual scale
    center_x = joint_cm(1) + joint_cm(3) / 2;                         % vertical center of first axis
    new_x0 = center_x - new_width_cm / 2;  

    marg_fh(end).Units = 'Centimeters';
    marg_fh(end).Position = [new_x0, joint_cm(2)+joint_cm(end)-0.5-4.5*(k==2), new_width_cm, 1, ];
    marg_pos = get(marg_fh(end),'Position'); % (5/5.2) = Corrected for the distortion caused by the color frame-circle
    set(marg_fh(end),'Position',marg_pos+[0 0.5 0 0])
    if k == 1
        set(marg_fh(end),'xaxislocation','top','XTickLabelRotation',0);

    end
    x_lbl = xlabel('$\dot{\theta}$ $[degrees/s]$  $\times$ $10^3$','Interpreter','latex');
    y_lbl = ylabel('Prob.');
    set(x_lbl,'Units','centimeters');
    if k == 2
        marg_fh(end).YDir = 'reverse';
        set(y_lbl,'rotation',270);
        set(y_lbl,'Units','centimeters');
        set(y_lbl,'Position',get(y_lbl,'Position')+[-0.6 0 0 ]);
    end
    marg_fh(end).XTickLabel = arrayfun(@(x) sprintf('%.0f', x / 1000), xticks, 'UniformOutput', false);

end
max_y_lim = max([marg_fh(end-1:end).YLim]);
ylim(marg_fh(end-1:end),[0,max_y_lim]);

%% Cont and non-cont alpha-distributions

% Make sure you're working on untouched alpha-vals
[alpha_vals,r_vals] = cart2pol(d_ang,d_curv);

% Plot a colorbar near the x-axis, for reference of alpha-values
interval_size = 0.1;
line_clr =  'rb';
fh(5) = axes();
x_clrs = flipud(copper(181));
x_tick_bounds = linspace(-pi,pi,length(x_clrs));
x_tick_vals = x_tick_bounds+0.5*mode(diff(x_tick_bounds));
x_val_dists = min(abs([circ_dist(x_tick_vals,alpha_star_val)',...
    circ_dist(x_tick_vals,-pi+alpha_star_val)']),[],2);
x_val_idx = discretize(x_val_dists,linspace(0,pi/2,length(x_clrs)));
x_clrs = x_clrs(x_val_idx,:);
rect_x = x_tick_bounds(1:end-1)';
rect_x(:,2) = abs(mode(diff(x_tick_bounds,1)));
rect_y = repmat([-0.0075 0.0075],size(rect_x,1),1);
rect_pos = [rect_x,rect_y];
rect_pos = rect_pos(:,[1,3,2,4]);
for k = 1:(length(x_tick_bounds)-1)
rectangle('Position',rect_pos(k,:),'FaceColor',x_clrs(k,:),'EdgeColor','None')
hold on
end

% Get and plot the alpha-distributions for contact and non-contact
dist_data = cell(0);
for k = 1:2
    if k == 1
        [y,x] = histcounts(alpha_vals(whObjDist<cont_thresh),-pi:interval_size:pi,'Normalization','pdf');
        dist_data{k} = alpha_vals(whObjDist<cont_thresh);
    else
        [y,x] = histcounts(alpha_vals(whObjDist>cont_thresh),-pi:interval_size:pi,'Normalization','pdf');
        dist_data{k} = alpha_vals(whObjDist>cont_thresh);
    end
    x = x(1:end-1)+0.5*mode(diff(x));
    y = y./sum(y);
    
    hold on
    plot([1 1]*pi/2,[0 0.1],'k')
    plot([1 1]*-pi/2,[0 0.1],'k')
    plot(x,y,'Color',line_clr(k));
end


% Format and beautify
fh(5).Units = 'Centimeters';
fh(5).Position = [15.5 23.5 3 1.25];
ax = axis();
ax(1:2) = [-pi pi];
axis(ax)
x_ticks = -pi:(pi/2):pi;
x_tick_labels = makePolarTicks(x_ticks);
set(gca,'XTick',x_ticks,'XTickLabels',x_tick_labels,...
    'XTickLabelRotation',45,'TickLabelInterpreter','Latex');
xlabel('$\alpha$ [radians]');
ylabel('Prob.');
x_label_handle = get(gca,'xlabel');
set(x_label_handle,'Units','centimeters','Interpreter','latex');
x_label_dist = get(x_label_handle,'position');
set(x_label_handle,'Position',x_label_dist+[0,-0.15,0]);
y_lim = ylim;
y_lim(2) = 0.06;
ylim(y_lim);

% Entropy permutation:
[alpha_vals,r_vals] = cart2pol(d_ang,d_curv);
N1 = sum(whObjDist<=cont_thresh);
N2 = sum(whObjDist>cont_thresh);
x_spacing = 361;
sim_res = [];
if get_permutation_p_val_flag
rng(121987,'twister') % Set random seed for reproducibility
for mm = 1:1001
    if mm ~= 1001
        shuff_idx = randperm(length(alpha_vals));
        grp1 = alpha_vals(shuff_idx(1:N1));
        grp2 = alpha_vals(shuff_idx(N1+1:end));
    else
        grp1 = alpha_vals(whObjDist<=cont_thresh);
        grp2 = alpha_vals(whObjDist>cont_thresh);
    end
    y1 = histcounts(grp1,linspace(-pi,pi,x_spacing),'Normalization','pdf');
    y2 = histcounts(grp2,linspace(-pi,pi,x_spacing),'Normalization','pdf');
    y1 = y1./sum(y1);
    y2 = y2./sum(y2);
    e1 = -nansum(y1.*log(y1));
    e2 = -nansum(y2.*log(y2));
    sim_res(mm) = (e2-e1);
end
ent_bstrp = mean(sim_res(1:end-1)>sim_res(end));
end

%% One dimensional heatmap - polar bar colored by distance (height by frequency)

% Make sure we're working on untouched alpha values
[alpha_vals,~] = cart2pol(d_ang,d_curv);

% Parameters:
max_dist_val = 0.16;
num_dist_clrs = 20;
num_alpha_categories = 361;
alpha_bounds = linspace(-pi,pi,num_alpha_categories);
dist_bounds = linspace(0,max_dist_val,num_dist_clrs);
alpha_categories = discretize(alpha_vals,alpha_bounds);
clrs = (jet(num_dist_clrs));

% Data collection
dists_for_heatmap = NaN(1,num_alpha_categories);
prob_of_alpha_val = NaN(1,num_alpha_categories);
for k = 1:(length(alpha_bounds)-1)
    dist_in_slice = whObjDist(alpha_categories == k);
    if length(dist_in_slice)<30
        continue
    end
    dists_for_heatmap(k) = sum(dist_in_slice<cont_thresh)./length(dist_in_slice);
    prob_of_alpha_val(k) = sum(alpha_categories == k)./length(alpha_categories);
end

% Plotting
fh(6) = axes();
for k = 1:(length(alpha_bounds)-1)
    x = alpha_bounds([k,k+1,k+1,k,k]);
    y = [0 0 1 1 0]*prob_of_alpha_val(k);
    [x,y] = pol2cart(x,y);
    if isnan(dists_for_heatmap(k))
        continue
    else
        [~,clr_idx] = min(abs(dists_for_heatmap(k)-dist_bounds));
        patch(x,y,clrs(clr_idx,:),'EdgeColor','None');
    end
end

circ_angles = linspace(-pi,pi,181);
circ_angle_for_clrs = linspace(0,pi/2,200);
circle_colormap = flipud(copper(200));
for i = 1:(length(circ_angles)-1)
    circ_angle_bounds = circ_angles(i:i+1);
    circ_angle_val = mean(circ_angle_bounds);
    circ_angle_val = min(...
                        abs(circ_dist(circ_angle_val,alpha_star_val)),...
                        abs(circ_dist(circ_angle_val,alpha_star_val-pi)));
    circ_angle_bounds = circ_angle_bounds([1,1,2,2,1]);
    circ_r_bounds = [0.006 0.0065];
    circ_r_bounds = circ_r_bounds([1 2 2 1 1]);
    [x,y] = pol2cart(circ_angle_bounds,circ_r_bounds);
    [~,idx] = min(abs(circ_angle_val-circ_angle_for_clrs));
    circ_clr = circle_colormap(idx,:);
    patch(x,y,circ_clr,'EdgeColor','None')
end

axis equal
axis(0.0065*[-1 1 -1 1])
mockPolarAxes(8,linspace(0,0.0065,5),fh(6),1.2,[0 0 0],0,{'0';'2';'4';'6';'8x10^{-3}'});
fh(6).Units = 'Centimeters';
fh(6).Position = [15 17 4.25 4.25];
colormap(fh(6),(jet(num_dist_clrs)));
clr_bar = colorbar(fh(6));
clr_bar.Units = 'Centimeters';
clr_bar.Position = [14.8 17.5 0.2 3.25];
practical_max_val = round(max(dists_for_heatmap),2);
practical_max_val_position = practical_max_val/max_dist_val;
set(clr_bar,'Ticks',[0,practical_max_val_position],'TickLabels',{'0';num2str(practical_max_val)});
clr_bar.Label.String = {'Contact Prob.'};
clr_bar.Label.Position = clr_bar.Label.Position+[2.6 0 0];
set(clr_bar,'TickDirection','Out')


%%  Probability of contact as a function of alpha
% Parameters:
alpha_categories_N = 11;
alpha_bounds = linspace(0,pi/2,alpha_categories_N);
alpha_categories_clipped = discretize(abs(alpha_vals-(alpha_star_val)),alpha_bounds);


% Make sure we're working on untouched alpha values
[alpha_vals,~] = cart2pol(d_ang,d_curv);

fh(7) = axes();

% Plot a colorbar near the x-axis, for reference of alpha-values
x_clrs = copper(181);
x_tick_bounds = linspace(pi,-pi,length(x_clrs));
x_tick_vals = x_tick_bounds+0.5*mode(diff(x_tick_bounds));
x_val_dists = min(abs([circ_dist(x_tick_vals,alpha_star_val)',...
circ_dist(x_tick_vals,-pi+alpha_star_val)']),[],2);
x_val_idx = discretize(x_val_dists,linspace(0,pi/2,length(x_clrs)));
x_clrs = x_clrs(x_val_idx,:);
rect_x = x_tick_bounds(1:end-1)';
rect_x(:,2) = abs(mode(diff(x_tick_bounds,1)));
rect_y = repmat([-0.02 0.02],size(rect_x,1),1);
rect_pos = [rect_x,rect_y];
rect_pos = rect_pos(:,[1,3,2,4]);
for k = 1:(length(x_tick_bounds)-1)
rectangle('Position',rect_pos(k,:),'FaceColor',x_clrs(k,:),'EdgeColor','None')
hold on
end


% Collect contact probability data and plot
all_xy = [];
scatter_clrs = jet(100);
clr_scale = linspace(0,max_dist_val,100);
for k = 1:max(alpha_categories_clipped)
    [y,x] = histcounts(whObjDist(alpha_categories_clipped == k),0:0.01:10,'Normalization','cdf');
    x = x(1:end-1)+0.5*mode(diff(x));
    [~,min_idx] = min(abs(x-cont_thresh));
    [~,clr_idx] = min(abs(y(min_idx)-clr_scale));
    sc = scatter(mean(alpha_bounds(k:k+1)'),y(min_idx),'o','filled','MarkerEdgeColor','k','MarkerFaceColor',scatter_clrs(clr_idx,:));
    hold on
    all_xy = [all_xy; [mean(alpha_bounds(k:k+1)'),y(min_idx)]];
end

% Fit and plot
% Initial guess from fit:
[a,~] = fit(all_xy(:,1),all_xy(:,2),'a*exp(-b*x)+c','StartPoint',[0.1 3 0.3]);
% Feed as initial guess into least-squares fit, to facilitate calculating
% statistics:
x0 = [a.a,a.b,a.c];
model = @(p, x) p(1) * exp(p(2) * x) + p(3);
x0 = [1, -1, 0];  % Initial guess
[beta_fit, ~, ~, ~, ~, ~, ~] = ...
    lsqcurvefit(model, x0, all_xy(:,1), all_xy(:,2));

% Get adjusted r-squared:
n = numel(all_xy(:,2));     % number of data points
k = numel(beta_fit);       % number of parameters
y_pred = model(beta_fit, all_xy(:,1));      % predicted values
residuals = all_xy(:,2) - y_pred;
SS_res = sum(residuals.^2);                      % Residual sum of squares
SS_tot = sum((all_xy(:,2) - mean(all_xy(:,2))).^2);        % Total sum of squares
R2 = 1 - SS_res / SS_tot;                        % R-squared
exp_fit_adj_r = 1 - (1 - R2)*(n - 1)/(n - k);            % Adjusted R-squared

% Get p-value vs constant
df1 = k - 1;
df2 = n - k;
MS_model = (SS_tot - SS_res) / df1;
MS_resid = SS_res / df2;
F = MS_model / MS_resid;
exp_fit_p_value = 1 - fcdf(F, df1, df2);

plot(all_xy(:,1),a(all_xy(:,1)),'k--','LineWidth',1);

% Format and beautify
fh(7).Units = 'Centimeters';
fh(7).Position = [15.5 14.45 3 1.25];
ylabel('Contact Prob.');
xlabel('$ \left|dist(\alpha^{*})\right| / {\frac{\pi}{2}}$','Interpreter','Latex');
x_ticks = linspace(0,pi/2,5);
set(gca,'XTickLabelRotation',0,'TickLabelInterpreter','Latex');
set(gca,'YTick',linspace(0,0.15,4));
set(gca, 'TickDir', 'out');

% Y-axis colormap
x_clrs = jet(20);
rect_y = linspace(0,max_dist_val,21)';
rect_y = [rect_y(1:end-1),diff(rect_y,1)];
rect_x = repmat([-0.075 0.075],size(rect_y,1),1);
rect_pos = [rect_x,rect_y];
rect_pos = rect_pos(:,[1,3,2,4]);
for k = 1:(size(rect_pos,1))
rectangle('Position',rect_pos(k,:),'FaceColor',x_clrs(k,:),'EdgeColor','None')
hold on
end
xlim([-0.075 pi/2])
ylim([-0.02 max_dist_val])
set(gca,'XTick',linspace(0,pi/2,5),'XTickLabel',(linspace(0,1,5)))

%% Tidy up
child_list = get(gcf,'children');
for k = 1:length(child_list)
    try
        set(child_list(k),'FontSize',10,'FontWeight','Normal');
    end
end
txt_list = findall(gcf,'Type','Text');
for k = 1:length(txt_list)
    set(txt_list(k),'FontSize',10,'FontWeight','Normal');
end
axes_list = findall(gcf,'Type','Axes');
for k = 1:length(axes_list)
    set(axes_list(k),'FontSize',10,'FontWeight','Normal');
    set(get(axes_list(k),'XLabel'),'FontSize',10,'FontWeight','Normal');
    set(get(axes_list(k),'YLabel'),'FontSize',10,'FontWeight','Normal');
end
clrbar_list = findall(gcf,'Type','Colorbar');
for k = 1:length(clrbar_list)
    set(clrbar_list(k),'FontSize',10,'FontWeight','Normal');
end
annot_list = findall(gcf,'Type','textboxshape');
for k = 1:length(annot_list)
    set(annot_list(k),'FontSize',10,'FontWeight','Normal');
end
drawnow expose

set([txt1,txt2,txt3],'FontSize',12)

%% Alternative phase-planes (extended figure)
% Do the same as Fig.2b-c above, separately for other pairs of variables
makeFullPagePDF();
supp_joint_fh = gobjects(0);
supp_marg_fh_ang = gobjects(0);
supp_marg_fh_curv = gobjects(0);
for var_iter = 1:3 % There's an option here to re-generate the original, but I skip it
    % It has slightly different values because the binnig strategy is
    % different, but the result is essentially the same

switch var_iter % Variable pairs
    case 1
        curv = cell2mat(cellfun(@(x) x.Tracks.smoothCurv(:)',allJobs,'Uni',0));
        curv_marginal_plot_range = linspace(-70,70,41);
        ang = cell2mat(cellfun(@(x) x.Tracks.smoothAngles(:)',allJobs,'Uni',0));
        ang_marginal_plot_range = linspace(-70,70,41);
        mot_lbl = '$|\theta|-90$ [${degrees}$]';
        sen_lbl = '${\kappa}$ [$m^{-1}$]';
        marg_ang_bounds = [0 0.06];
        marg_curv_bounds = [0 0.08];
    case 2
        curv = cell2mat(cellfun(@(x) (500/1000)*x.Tracks.diffCurv(:)',allJobs,'Uni',0));
        curv_marginal_plot_range = linspace(-1.5,1.5,41);
        ang = cell2mat(cellfun(@(x) x.Tracks.smoothAngles(:)',allJobs,'Uni',0));
        ang_marginal_plot_range = linspace(-70,70,41);
        mot_lbl = '$|\theta|-90|$ [${degrees}$]';
        sen_lbl = '$\dot{\kappa}$ $[m^{-1}/s]$ $\times$ $10^3$';
        marg_ang_bounds = [0 0.06];
        marg_curv_bounds = [0 0.1];
    case 3
        curv = cell2mat(cellfun(@(x) x.Tracks.smoothCurv(:)',allJobs,'Uni',0));
        curv_marginal_plot_range = linspace(-70,70,41);
        ang = cell2mat(cellfun(@(x) (500/1000)*x.Tracks.diffAngles(:)',allJobs,'Uni',0));
        ang_marginal_plot_range = linspace(-1.5,1.5,41);
        mot_lbl = '$\dot{\theta}$ $[degrees/s]$  $\times$ $10^3$';
        sen_lbl = '${\kappa}$ [$m^{-1}$]';
        marg_ang_bounds = [0 0.14];
        marg_curv_bounds = [0 0.08];
    case 4
        curv = cell2mat(cellfun(@(x) (500/1000)*x.Tracks.diffCurv(:)',allJobs,'Uni',0));
        curv_marginal_plot_range = linspace(-1.5,1.5,41);
        ang = cell2mat(cellfun(@(x) (500/1000)*x.Tracks.diffAngles(:)',allJobs,'Uni',0));
        ang_marginal_plot_range = linspace(-1.5,1.5,41);
        mot_lbl = '$\dot{\theta}$ $[degrees/s]$  $\times$ $10^3$';
        sen_lbl = '$\dot{\kappa}$ $[m^{-1}/s]$ $\times$ $10^3$';
        marg_ang_bounds = [0 0.14];
        marg_curv_bounds = [0 0.1];
end

% No snout-object contact frames
rawWhObjDist = cell2mat(cellfun(@(x) x.Tracks.whiskerObjDist(:)',allJobs,'Uni',0));
rawWhObjDist = rawWhObjDist./(px2mm*10);
snObjDist = cell2mat(cellfun(@(x) x.Tracks.snoutObjDist(:)',allJobs,'Uni',0));
snObjDist = cellfun(@(x) x.Tracks.snoutObjDist,allJobs,'Uni',0);
n_row = cellfun(@(x) size(x.Tracks.whiskerObjDist,1),allJobs,'Uni',0);
snObjDist = cellfun(@(x,y) repmat(x,y,1),snObjDist,n_row,'Uni',0);
snObjDist = cell2mat(cellfun(@(x) x(:)',snObjDist,'Uni',0));
snObjDist = snObjDist./(px2mm*10);
if excludeNoseTipCont
    rawWhObjDist(snObjDist<noseTipContThresh) = NaN;
end
pre_cont_dist = cellfun(@(x) repmat((cumsum(any((x.Tracks.whiskerObjDist./(px2mm*10)) < pre_cont_thresh, 1), 2) == 0),size(x.Tracks.whiskerObjDist,1),1), allJobs, 'Uni',0);
pre_cont_dist = cell2mat(cellfun(@(x) x(:)',pre_cont_dist,'Uni',0));

% For raw theta, center the distribution around 90deg
if var_iter < 3 
    ang = abs(ang)-90;
end

% Get Ns of observations for summary statistics
nan_vals = isnan(rawWhObjDist) | isnan(ang) | isnan(curv);
curv = curv(~nan_vals);
ang = ang(~nan_vals);
whObjDist = rawWhObjDist(~nan_vals);
pre_cont_dist = pre_cont_dist(~nan_vals);
N_pre = sum(pre_cont_dist);
N_non = sum(~(pre_cont_dist | (whObjDist<=cont_thresh)));
N_cont = sum(~(pre_cont_dist | (whObjDist>cont_thresh)));
Ns_raw_and_dot(var_iter,:) = [N_pre,N_cont,N_non];

% Making it easier to plot joint when units are very different:
max_mot_var = max(ang_marginal_plot_range);
max_sen_var = max(curv_marginal_plot_range);
ang = ang./max_mot_var;
curv = curv./max_sen_var;
max_r = 1;

% 2D visit-rates: All, cont and non-cont
r_joint_plot_range = linspace(0,max_r,length(r_joint_plot_range));
% Gather count data
counts_for_heatmap = cell(1,3);
for mm = 1:3
    [alpha_vals,r_vals] = cart2pol(ang,curv);
    
if mm == 1
    if pre_cont_flag
        alpha_vals(~pre_cont_dist) = [];
        r_vals(~pre_cont_dist) = [];
    end
elseif mm == 2
    alpha_vals(pre_cont_dist | (whObjDist>cont_thresh)) = [];
    r_vals(pre_cont_dist | (whObjDist>cont_thresh)) = [];
elseif mm == 3
    alpha_vals(pre_cont_dist | (whObjDist<=cont_thresh)) = [];
    r_vals(pre_cont_dist | (whObjDist<=cont_thresh)) = [];
end
alpha_vals_binned = discretize(alpha_vals,alpha_joint_plot_range);
r_val_binned = discretize(r_vals,r_joint_plot_range);

counts_for_heatmap{mm} = NaN(n_polar_wedges-2,n_r_circles-2);
for k = 1:(length(alpha_joint_plot_range)-1)
    for kk = 1:(length(r_joint_plot_range)-1)
        counts_for_heatmap{mm}(k,kk) = sum(alpha_vals_binned == k & r_val_binned == kk);
    end
end
end

% Get probabilities from counts:
probs_for_heatmap = cellfun(@(x) x./nansum(x(:)),counts_for_heatmap,'Uni',0);
log_scale_clrs = true;
prob_bins = logspace(-4,-2.75,21);
prob_clrs = gray(length(prob_bins)+1);
for mm = 1:3
supp_joint_fh(var_iter,mm) = axes();
for k = 1:(length(alpha_joint_plot_range)-1)
    for kk = 1:(length(r_joint_plot_range)-1)

        x = alpha_joint_plot_range([k,k+1,k+1,k,k]);
        y = r_joint_plot_range([kk,kk,kk+1,kk+1,kk]);
        [x,y] = pol2cart(x,y);
        if isnan(probs_for_heatmap{mm}(k,kk))
            clr_idx = 1;
        else
            if log_scale_clrs
                [~,clr_idx] = min(abs((probs_for_heatmap{mm}(k,kk))-prob_bins));
            end
        end
        patch(x,y,prob_clrs(clr_idx,:),'EdgeColor','None');
    end
end

        
% Add colored borders
max_r_val_diff = max(r_joint_plot_range);
border_circ_alpha = linspace(-pi,pi,37);
border_circ_alpha = [border_circ_alpha(1:end-1);border_circ_alpha(2:end)]';
border_circ_r = (max_r_val_diff.*[1 1]+[0 0.075]).*ones(size(border_circ_alpha));
if mm == 1
    circ_clr = [0 0 1];
elseif mm == 2
    circ_clr = [1 0 0];
elseif mm == 3
    circ_clr = [0 0 1];
end
border_clrs = repmat(circ_clr,size(border_circ_alpha,1),1);
polarHeatMap(border_circ_alpha,border_circ_r,border_clrs);


% Adjust axes, ticks, labels
supp_joint_fh(var_iter,mm).Units = 'Centimeters';
set(supp_joint_fh(var_iter,mm),'TickLength',[0 0],'xticklabels',{},'yticklabels',{});
set(gca,'Color','None','XColor','None','YColor','None');
axis([-1 1 -1 1]*(max_r_val_diff+0.2))
axis equal
if (var_iter) == 1 & (mm==1)
    
    colormap(prob_clrs);
    clr_bar = colorbar('Southoutside');
    clr_bin_centers = linspace(0,1,length(prob_clrs));
    clr_bin_centers = clr_bin_centers(1:end-1)+0.5*mode(diff(clr_bin_centers));
    if log_scale_clrs
        y_ticks = interp1(prob_bins, clr_bin_centers, 10.^(clr_ticks), 'linear', 'extrap');
        y_lbls =  arrayfun(@(x) ['10^{',num2str(x),'}'],clr_ticks,'Uni',0);
    end
    set(clr_bar,'YTick',y_ticks,'YTickLabel',y_lbls);
    set(clr_bar,'Units','centimeters','Position',[1.75 21 3.5 0.2]);
    clr_bar.Label.String = 'Probability Mass';

end


% Adjust positions

max_r_val_diff = max(r_joint_plot_range);
if mm == 1
    supp_joint_fh(var_iter,mm).InnerPosition = [1.75 17-7*(var_iter-1) 3.25 3.25];
elseif mm == 2
    supp_joint_fh(var_iter,mm).InnerPosition = [7.5 17-7*(var_iter-1) 3.25 3.25];
elseif mm == 3
    supp_joint_fh(var_iter,mm).InnerPosition = [13.25 17-7*(var_iter-1) 3.25 3.25];
end
axes(supp_joint_fh(var_iter,mm))
axis([-1 1 -1 1]*(max_r_val_diff+0.075))

end

if var_iter == 1
txt1 = text(supp_joint_fh(var_iter,1), 0.5, 1.8, 'Pre-contact', 'Units', 'normalized', ...
     'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', 'FontSize', 12);
txt2 = text(supp_joint_fh(var_iter,2), 0.5, 1.8, 'Contact', 'Units', 'normalized', ...
     'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', 'FontSize', 12);
txt3 = text(supp_joint_fh(var_iter,3), 0.5, 1.7, {'Non-contact';'(after 1st cont.)'}, 'Units', 'normalized', ...
     'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', 'FontSize', 12);
end
axes(supp_joint_fh(var_iter,3))
% Going back original units for marginals
ang = ang.*max_mot_var;
curv = curv.*max_sen_var;
% Motor variable marginals:
for mm = 2:3
    if mm == 2
        marg_ang = ang(~(pre_cont_dist | (whObjDist>cont_thresh)));
        marg_curv = curv(~(pre_cont_dist | (whObjDist>cont_thresh)));
    elseif mm == 3
        marg_ang = ang(~(pre_cont_dist | (whObjDist<=cont_thresh)));
        marg_curv = curv(~(pre_cont_dist | (whObjDist<=cont_thresh)));
    end
    
    supp_marg_fh_ang(var_iter,mm-1) = axes('Units','centimeters');
    [y,x] = histcounts(marg_ang,ang_marginal_plot_range,'Normalization','pdf');
    y = y./sum(y(:));
    x = x(1:end-1)+0.5*mode(diff(x));
    if mm == 2
        plot(x,y,'Color','r');
    elseif mm == 3
        plot(x,y,'Color','b');
    end
    hold on
    xlim([min(ang_marginal_plot_range),max(ang_marginal_plot_range)])
    joint_x_lim = get(supp_joint_fh(var_iter,mm),'XLim');
    joint_cm =  get(supp_joint_fh(var_iter,mm),'Position');
    scale = (2 * joint_x_lim(2)) / joint_cm(3);                  % data units per cm
    new_width_cm = (2 *  max(r_joint_plot_range)) / scale;           % height needed for same visual scale
    center_x = joint_cm(1) + joint_cm(3) / 2;                         % vertical center of first axis
    new_x0 = center_x - new_width_cm / 2;  
    supp_marg_fh_ang(var_iter,mm-1).Position = [new_x0, joint_cm(2)+joint_cm(end)-0.5, new_width_cm, 1, ];
    marg_pos = get(supp_marg_fh_ang(var_iter,mm-1),'Position'); % (5/5.2) = Corrected for the distortion caused by the color frame-circle
    set(supp_marg_fh_ang(var_iter,mm-1),'Position',marg_pos+[0 0.5 0 0])
    set(supp_marg_fh_ang(var_iter,mm-1),'xaxislocation','top','XTickLabelRotation',0);
    x_lbl = xlabel(mot_lbl,'Interpreter','latex');
    y_lbl = ylabel('Prob.');
    set(x_lbl,'Units','centimeters');
end
ylim([supp_marg_fh_ang(var_iter,:)],marg_ang_bounds);
% Sensory variable marginals:
for mm = 2:3
    if mm == 2
        marg_ang = ang(~(pre_cont_dist | (whObjDist>cont_thresh)));
        marg_curv = curv(~(pre_cont_dist | (whObjDist>cont_thresh)));
    elseif mm == 3
        marg_ang = ang(~(pre_cont_dist | (whObjDist<=cont_thresh)));
        marg_curv = curv(~(pre_cont_dist | (whObjDist<=cont_thresh)));
    end


    supp_marg_fh_curv(var_iter,mm-1) = axes('Units','centimeters');
    [y,x] = histcounts(marg_curv,curv_marginal_plot_range,'Normalization','pdf');
    y = y./sum(y(:));
    x = x(1:end-1)+0.5*mode(diff(x));
    if mm == 2
        plot(y,x,'Color','r');
    elseif mm == 3
        plot(y,x,'Color','b');
    end
    hold on
    ylim([min(curv_marginal_plot_range),max(curv_marginal_plot_range)])
    joint_y_lim = get(supp_joint_fh(var_iter,mm),'YLim');
    joint_cm =  get(supp_joint_fh(var_iter,mm),'Position');
    scale = (2 * joint_y_lim(2)) / joint_cm(end);                  % data units per cm
    new_height_cm = (2 *  max(r_joint_plot_range)) / scale;           % height needed for same visual scale
    center_y = joint_cm(2) + joint_cm(end) / 2;                         % vertical center of first axis
    new_y0 = center_y - new_height_cm / 2;                  % align centers
    supp_marg_fh_curv(var_iter,mm-1).Units = 'Centimeters';
    supp_marg_fh_curv(var_iter,mm-1).Position = [10.75+(mm==3)*5.75, new_y0, 1, new_height_cm];
    xlabel('Prob.');
    set(supp_marg_fh_curv(var_iter,mm-1),'yaxislocation','right');
    set(supp_marg_fh_curv(var_iter,mm-1),'YTickLabelRotation',-45);
    y_lbl = ylabel(sen_lbl,'rotation',90,'Interpreter','latex');
    set(y_lbl,'Units','centimeters');
    set(y_lbl,'Position',get(y_lbl,'Position')+[0 0 0 ]);
    y_ticks = yticks();

    y_ticks = get(gca,'YTick');
    new_y_tick_labels = y_ticks./curv_scale_coeff;
    scale_log_ten = floor(log10(max(abs(new_y_tick_labels))));
    new_y_tick_labels = round(new_y_tick_labels./(10^scale_log_ten));
    new_y_tick_labels = new_y_tick_labels.*(10^scale_log_ten);
    new_y_ticks = new_y_tick_labels.*curv_scale_coeff;
    set(gca,'YTick',new_y_ticks,'YTickLabel',new_y_tick_labels);
    

end
xlim([supp_marg_fh_curv(var_iter,:)],marg_curv_bounds);

% Late addition, add axis ticks for pre-contact:
axes(supp_joint_fh(var_iter,1));
set(gca,'Xcolor','k','YColor','k');
plotted_yticks = supp_marg_fh_curv(var_iter,1).YTick./max(supp_marg_fh_curv(var_iter,1).YLim);
set(gca,'YTick',plotted_yticks,'YTickLabel',supp_marg_fh_curv(var_iter,1).YTickLabel);
plotted_xticks = supp_marg_fh_ang(var_iter,1).XTick./max(supp_marg_fh_ang(var_iter,1).XLim);
set(gca,'XTick',plotted_xticks,'XTickLabel',supp_marg_fh_ang(var_iter,1).XTickLabel);
set(gca,'TickLength',supp_marg_fh_curv(1).TickLength)
xlabel(mot_lbl,'Interpreter','latex');
ylabel(sen_lbl,'Interpreter','latex');

% Tidy up each figure:
child_list = get(gcf,'children');
for k = 1:length(child_list)
    try
        set(child_list(k),'FontSize',10,'FontWeight','Normal');
    end
end
txt_list = findall(gcf,'Type','Text');
for k = 1:length(txt_list)
    set(txt_list(k),'FontSize',10,'FontWeight','Normal');
end
axes_list = findall(gcf,'Type','Axes');
for k = 1:length(axes_list)
    set(axes_list(k),'FontSize',10,'FontWeight','Normal');
    set(get(axes_list(k),'XLabel'),'FontSize',10,'FontWeight','Normal');
    set(get(axes_list(k),'YLabel'),'FontSize',10,'FontWeight','Normal');
end
clrbar_list = findall(gcf,'Type','Colorbar');
for k = 1:length(clrbar_list)
    set(clrbar_list(k),'FontSize',10,'FontWeight','Normal');
end
annot_list = findall(gcf,'Type','textboxshape');
for k = 1:length(annot_list)
    set(annot_list(k),'FontSize',10,'FontWeight','Normal');
end
set([txt1,txt2,txt3],'FontSize',12)
drawnow expose
end


%% Report statistics:
disp('Fig 2 stats:')
disp('Fig 2a-c: Ns dot-variables, pre-cont, cont and non-cont:')
disp(main_fig_overall_Ns)

disp('extended figure: Ns raw-variables, pre-cont, cont and non-cont:')
helper_ttls = {'Both variables raw';'Theta, Kappa-dot';'Theta-dot, Kappa'};
for k = 1:size(Ns_raw_and_dot,1)
    disp([helper_ttls{k},':']);
    disp(Ns_raw_and_dot(k,:));
end


if get_permutation_p_val_flag
disp(' ')
disp('Fig2a-c: p(JSD(cont,non-cont)-JSD(non-cont,pre-cont) <= label-shuffled null)')
disp('Tested with n-iterations = 1000')
disp(num2str(p_val_JSD_diffs))

disp(' ')
disp('Fig 2d: p(entropy(cont)-entropy(non-cont) <= label-shuffled null)')
disp('Tested with n-iterations = 1000')
disp(num2str(ent_bstrp))

end

disp(' ')
disp('Fig 2f: p(cont) as function of a^* distance:')
disp(['adj. r^2 = ',num2str(exp_fit_adj_r),'; p = ',num2str(exp_fit_p_value)])

end