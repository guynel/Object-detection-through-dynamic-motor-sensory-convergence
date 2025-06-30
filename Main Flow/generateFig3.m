function [allJobs] = generateFig3(allJobs,px2mm,step_N, max_step,extend_to_cycle_end_flag,cont_thresh_px,alpha_star_displacement)
get_permutation_p_val_flag = 1; % Do you want to do permutation test? Resoure-heavy...
%% Gather contact and non-contact flow data

% Pre-allocate
bin_num = 121;
all_cont_alpha = [];
all_non_cont_alpha = [];
all_cont_d_alpha = [];
all_non_cont_d_alpha = [];

% Collect data
for trial_num = 1:length(allJobs)
    trim_idx = allJobs{trial_num}.Tracks.trim_idx;
    if isnan(trim_idx)
        continue
    end
    dd = allJobs{trial_num}.Tracks.whiskerObjDist(:,1:trim_idx)<cont_thresh_px;
    cont = dd./dd;
    non_cont = (~dd)./(~dd);

    d_ang = allJobs{trial_num}.Tracks.diffAngles(:,1:trim_idx);
    d_curv = allJobs{trial_num}.Tracks.diffCurv(:,1:trim_idx);
    [alpha_vals,r_vals] = cart2pol(d_ang,d_curv);
    cont_alpha = alpha_vals.*cont;
    non_cont_alpha = alpha_vals.*non_cont;
    % Just general statistics within approach only:
    all_cont_alpha = [all_cont_alpha; cont_alpha(~isnan(cont_alpha))];
    all_non_cont_alpha = [all_non_cont_alpha; non_cont_alpha(~isnan(non_cont_alpha))];
    % Timeseries pairs for cont:
    d_alpha = circ_dist(cont_alpha(:,step_N:end),cont_alpha(:,1:(end-(step_N-1))));
    origin_alpha = cont_alpha(:,1:(end-(step_N-1)));
    d_alpha = [origin_alpha(:),d_alpha(:)];
    d_alpha = d_alpha(~any(isnan(d_alpha),2),:);
    all_cont_d_alpha = [all_cont_d_alpha; d_alpha];

    % Timeseries pairs for non-cont:
    d_alpha = circ_dist(non_cont_alpha(:,step_N:end),non_cont_alpha(:,1:(end-(step_N-1))));
    origin_alpha = non_cont_alpha(:,1:(end-(step_N-1)));
    d_alpha = [origin_alpha(:),d_alpha(:)];
    d_alpha = d_alpha(~any(isnan(d_alpha),2),:);
    all_non_cont_d_alpha = [all_non_cont_d_alpha; d_alpha];
    
end
flow_data = {all_cont_d_alpha;all_non_cont_d_alpha};

%% Start with control:
% Resample from data of the same size as contact, from the entire population
% i.e. label-collapsed data; repeat 100 times; will be used as control
rng(121987,'twister')
for iter = 1:100
    pooled_datapts = [flow_data{1}; flow_data{2}];
    pooled_datapts(abs(pooled_datapts(:,2))>max_step,:) = [];
    control_alpha_idx = randsample(1:length(pooled_datapts),length(flow_data{1}),1);
    control_alpha = pooled_datapts(control_alpha_idx,:);

    h = 0.15;

    % Plot


    alpha_bounds = linspace(-pi,pi,bin_num);
    alpha_centers = alpha_bounds(1:end-1)+0.5*mode(diff(alpha_bounds));

    data_vec = control_alpha;
    
    % Get directional tendency function, after kernel-smoothing:
    krnl = @(xi, x0, h) exp(-0.5 * ((xi - x0) ./ h).^2) ./ (h * sqrt(2 * pi));
    x0 = alpha_centers;
    kernel_func = zeros(size(x0));
    for i = 1:length(x0)
        weights = krnl(data_vec(:,1), x0(i), h);  % Compute weights for all x_i
        weights = weights./sum(weights);
        kernel_func(i) = -0.5+sum((data_vec(:,2) > 0) .* weights); % Weighted average of Y
    end
    smoothed_direction_control(iter,:) = kernel_func;
    
    % Get mobility function, after kernel-smoothing:
    kernel_func = zeros(size(x0));
    for i = 1:length(x0)
        weights = krnl(data_vec(:,1), x0(i), h);  % Compute weights for all x_i
        weights = weights./sum(weights);
        abs_data = abs(data_vec(:,2));
        [sorted_data, idx] = sort(abs_data);
        sorted_weights = weights(idx);
        cumsum_weights = cumsum(sorted_weights);
        med_idx = find(cumsum_weights >= 0.5, 1);  % Find the median index
        kernel_func(i) = sorted_data(med_idx);  % The weighted median absolute value
    end
    smoothed_mobility_control(iter,:) = kernel_func;


end


%% Caclulate kernel on actual data
% Repeats the process from above
h = 0.15;
flow_data = {all_cont_d_alpha;all_non_cont_d_alpha};

for k = 1:length(flow_data)
    data_vec = flow_data{k};
    data_vec(abs(data_vec(:,2))>max_step,:) = [];
    flow_data{k} = data_vec;
end

smoothed_direction_vec = cell(0);
smoothed_mobility_vec = cell(0);
for k = 1:length(flow_data)
    data_vec = flow_data{k};
    krnl = @(xi, x0, h) exp(-0.5 * ((xi - x0) ./ h).^2) ./ (h * sqrt(2 * pi));
    x0 = alpha_centers;
    kernel_func = zeros(size(x0));

    for i = 1:length(x0)
        weights = krnl(data_vec(:,1), x0(i), h);  % Compute weights for all x_i
        weights = weights./sum(weights);
        kernel_func(i) = -0.5+sum((data_vec(:,2) > 0) .* weights); % Weighted average of Y
    end
    smoothed_direction_vec{k} = kernel_func;

    kernel_func = zeros(size(x0));
    for i = 1:length(x0)
        weights = krnl(data_vec(:,1), x0(i), h);  % Compute weights for all x_i
        weights = weights./sum(weights);
        abs_data = abs(data_vec(:,2));
        [sorted_data, idx] = sort(abs_data);
        sorted_weights = weights(idx);
        cumsum_weights = cumsum(sorted_weights);
        med_idx = find(cumsum_weights >= 0.5, 1);  % Find the median index
        kernel_func(i) = sorted_data(med_idx);  % The weighted median absolute value
    end
    smoothed_mobility_vec{k} = kernel_func;
end
drawnow

% 
smoothed_mobility_vec = [smoothed_mobility_vec, mean(smoothed_mobility_control)];
smoothed_direction_vec = [smoothed_direction_vec, mean(smoothed_direction_control)];

%% Gather errorbars data
% Same process, re-sampling each group with replacement, repeated 100 times
flow_data = {all_cont_d_alpha;all_non_cont_d_alpha};
for k = 1:length(flow_data)
    data_vec = flow_data{k};
    data_vec(abs(data_vec(:,2))>max_step,:) = [];
    flow_data{k} = data_vec;
end

h = 0.15;

smoothed_mob_bstrp = cell(0);
smoothed_dir_bstrp = cell(0);
for k = 1:length(flow_data)
    bstrp_dir = [];
    bstrp_mob = [];
    for iter = 1:100
        data_vec = flow_data{k};
        data_vec(abs(data_vec(:,2))>max_step,:) = [];
        data_vec = data_vec(randsample(length(data_vec),length(data_vec),1),:);
        krnl = @(xi, x0, h) exp(-0.5 * ((xi - x0) ./ h).^2) ./ (h * sqrt(2 * pi));
        x0 = alpha_centers;
        kernel_func = zeros(size(x0));

        for i = 1:length(x0)
            weights = krnl(data_vec(:,1), x0(i), h);  % Compute weights for all x_i
            weights = weights./sum(weights);
            kernel_func(i) = -0.5+sum((data_vec(:,2) > 0) .* weights); % Weighted average of Y
        end
        bstrp_dir(iter,:) = kernel_func;

        kernel_func = zeros(size(x0));
        for i = 1:length(x0)
            weights = krnl(data_vec(:,1), x0(i), h);  % Compute weights for all x_i
            weights = weights./sum(weights);
            abs_data = abs(data_vec(:,2));
            [sorted_data, idx] = sort(abs_data);
            sorted_weights = weights(idx);
            cumsum_weights = cumsum(sorted_weights);
            med_idx = find(cumsum_weights >= 0.5, 1);  % Find the median index
            kernel_func(i) = sorted_data(med_idx);  % The weighted median absolute value
        end
        bstrp_mob(iter,:) = kernel_func;
    end
    % Take percentiles:
    smoothed_mob_bstrp{k} = prctile(bstrp_mob,[2.5,97.5]);
    smoothed_dir_bstrp{k} = prctile(bstrp_dir,[2.5,97.5]);
end
smoothed_mob_bstrp{3} = prctile(smoothed_mobility_control,[2.5,97.5]);
smoothed_dir_bstrp{3} = prctile(smoothed_direction_control,[2.5,97.5]);

%% Begin by plotting both mobility and directional tendency, as arrows

clrs = [1 0 0; 0 0 1];
main_fig = makeFullPagePDF();
h = gobjects(0);
h(1) = axes('Units','Centimeters','Position',[2.5 24.5 9 1]);
plot(NaN,NaN,'r','LineWidth',3);
hold on
plot(NaN,NaN,'b','LineWidth',3);
for mm = 1:2
    clr = clrs(mm,:);
    
    mob_vec = 0.5*smoothed_mobility_vec{mm}; % Scaling doesn't matter here

    prob_directional_tendency = smoothed_direction_vec{mm};
    max_dir = max(cellfun(@(x) max(abs(x)),smoothed_direction_vec));
    prob_directional_tendency = prob_directional_tendency./max_dir;
    for k = 2:3:(length(prob_directional_tendency)-1)
        x = alpha_centers(k)+0.5*diff(alpha_centers(k:k+1));
        
        arrow_vec = [0.15 0; 0.15 0.5; 0.5 0.5; 0 1; -0.5 0.5; -0.15 0.5; -0.15 0];
        angle = -pi/2;
        R = [cos(angle), -sin(angle); sin(angle), cos(angle)];
        arrow_vec = (arrow_vec *  R');
        if prob_directional_tendency(k)<0
            arrow_vec(:,1) = -1*arrow_vec(:,1);
        end
        arrow_vec = [arrow_vec; arrow_vec(1,:)];

        scaled_arrow = arrow_vec*(0.15+abs(mob_vec(k)));
        scaled_arrow = 0.35.*[0.8 1].*arrow_vec*(0.15+abs(mob_vec(k)));
        scaled_arrow = scaled_arrow+[x,1+(mm-1)*0.2];
        fill(scaled_arrow(:,1), scaled_arrow(:,2), clr, 'EdgeColor', 'none','FaceAlpha',abs(prob_directional_tendency(k)));
        hold on

    end

    plot(-alpha_star_displacement+[1 1]*pi/2,ylim,'k--')
    plot(-alpha_star_displacement+[1 1]*-pi/2,ylim,'k--')
    xlim([-pi pi])
    xlim([-pi pi]);
    set(gca,'XTickLabel',[],'YTick',[]);
    ylabel('Flow');
end
set(gca,'box','off','XColor','k','YColor','k');


%% Next plot mobility function

clrs = [1 0 0; 0 0 1; 0.5 0.5 0.5; 0.5 0.5 0.5; 0.5 0.5 0.5];
h(2) = axes('Units','Centimeters','Position',[2.5 21.9 9 2]);
for mm = 1:3
    clr = clrs(mm,:);
    
    alpha_dot_median = smoothed_mobility_vec{mm};
    if mm ~= 3
        plot(alpha_centers,alpha_dot_median,'Color',clr,'LineWidth',1)
    end
    hold on
    xx = [alpha_centers,fliplr(alpha_centers)];
    yy = smoothed_mob_bstrp{mm};
    yy = [yy(1,:),fliplr(yy(2,:))];
    if mm~=2
        patch(xx,yy(:),clr,'FaceAlpha',0.4,'EdgeColor','None');
    end
    
   

    plot([1 1]*(-pi/2-alpha_star_displacement),[-1,1],'--','Color',[0 0 0])
    plot([1 1]*(pi/2-alpha_star_displacement),[-1,1],'--','Color',[0 0 0])
    xlim([-pi pi])
    xlim([-pi pi]);
    set(gca,'XTick',[]);
    ylabel({'Mobility';'$median(|\Delta{\alpha}|)$'},'Interpreter','latex');
end
set(gca,'YTick',linspace(0,pi/2,5),'YTickLabel',makePolarTicks(linspace(0,pi/2,5)),'TickLabelInterpreter','latex');
set(gca,'box','off','XColor','k','YColor','k');
ylim([0 0.9]);
xlim([-pi pi]);

%% Next plot directional tendency function

clrs = [1 0 0; 0 0 1; 0.5 0.5 0.5];
h(3) = axes('Units','Centimeters','Position',[2.5 19.2 9 2]);

hold on

for mm = 1:3
    clr = clrs(mm,:);
    
    if mm ~= 3
        plot(alpha_centers,smoothed_direction_vec{mm},'Color',clr,'LineWidth',1)
    end
    hold on
    xx = [alpha_centers,fliplr(alpha_centers)];
    yy = smoothed_dir_bstrp{mm};
    yy = [yy(1,:),fliplr(yy(2,:))];
    if mm~=2
        patch(xx,yy(:),clr,'FaceAlpha',0.4,'EdgeColor','None');
    end
   
    plot([1 1]*(-pi/2-alpha_star_displacement),[-1,1],'--','Color',[0 0 0])
    plot([1 1]*(pi/2-alpha_star_displacement),[-1,1],'--','Color',[0 0 0])
    xlim([-pi pi])
    xlim([-pi pi]);
    x_ticks = [-pi,-pi/2-alpha_star_displacement,0,pi/2-alpha_star_displacement,pi];
    x_tick_lbls = makePolarTicks(linspace(-pi,pi,3));
    x_tick_lbls = {x_tick_lbls{1},'$-\alpha^{*}$',x_tick_lbls{2},'$\alpha^{*}$',x_tick_lbls{3}};
    set(gca,'XTick',x_ticks,'XTickLabel',x_tick_lbls,'TickLabelInterpreter','latex');
    xlabel('\alpha [radians]');
    ylabel({'Directional';'Tendency';'$p(\Delta\alpha > 0)-0.5$'},'Interpreter','latex');
end
set(gca,'box','off','XColor','k','YColor','k');
set(gca,'YTick',[-1 1]*0.15)

ylim([-1 1]*0.2)
xlim([-pi pi])
plot(xlim,[0 0],'k--')

%% Next plot arrows again, but on circle
% Here we also use the control
h(4) = axes('Units','Centimeters','Position',[12.45 19.7 5.3 5.3]);
hold on
% Some workaround to get a nice legend
plot(NaN,NaN,'r','LineWidth',3);
plot(NaN,NaN,'b','LineWidth',3);
plot(NaN,NaN,'Color',[0.5 0.5 0.5],'LineWidth',3);
mockPolarAxes(12,[0 1.5],gca,1.2,[0 0 0],0);
% Fix the polar-ticks for nice presentation:
txt_list = findall(gca,'Type','Text');
txt_list(logical(strcmpi({txt_list.String},'1.5'))).String = '';
txt_list(logical(strcmpi({txt_list.String},'$\frac{1}{3}\pi$'))).String = '';

% Plot alpha-star line
hold on
[x,y] =  pol2cart([-1 1]*pi/2-alpha_star_displacement,[1,1]*1.5);
plot(x,y,'--','Color',[0 0 0],'LineWidth',2)

% Plot arrows:
for mm = 1:length(smoothed_mobility_vec)
    clr = clrs(mm,:);
    mob_vec = 0.5*smoothed_mobility_vec{mm};
    prob_directional_tendency = smoothed_direction_vec{mm};
    max_dir = max(cellfun(@(x) max(abs(x)),smoothed_direction_vec));
    prob_directional_tendency = prob_directional_tendency./max_dir;

    for k = 2:5:(length(prob_directional_tendency)-1)
        x = alpha_centers(k);
        y = 1.5-(0.3*(mm-1));
        if prob_directional_tendency(k)<0
            u = alpha_centers(k-1);
            v =  1.5-(0.3*(mm-1));
        elseif prob_directional_tendency(k)>0
            u = alpha_centers(k+1);
            v =  1.5-(0.3*(mm-1));
        end
        [x,y] = pol2cart(x,y);
        [u,v] = pol2cart(u,v);

        arrow_vec = [0.15 0; 0.15 0.5; 0.5 0.5; 0 1; -0.5 0.5; -0.15 0.5; -0.15 0];
        angle = -pi/2;
        R = [cos(angle), -sin(angle); sin(angle), cos(angle)];
        arrow_vec = (arrow_vec *  R');
        arrow_vec = [arrow_vec; arrow_vec(1,:)];

        dx = (u-x);
        dy = (v-y);
        angle = atan2(dy, dx);
        R = [cos(angle), -sin(angle); sin(angle), cos(angle)];
        scaled_arrow = 0.7*arrow_vec*(0.15+0.75*abs(mob_vec(k)));
        scaled_arrow = (scaled_arrow *  R')+[x,y];
      
        fill(scaled_arrow(:,1), scaled_arrow(:,2), clr, 'EdgeColor', 'none','FaceAlpha',abs(prob_directional_tendency(k)));
        hold on
        axis equal
        axis([-1.1 1.1 -1.1 1.1]*1.5)
        set(gca,'XTick',[],'YTick',[],'XTickLabel','','YTickLabel','','box','off','XColor','None','YColor','None');
    end

end

%% Add protraction/retraction bounds as illustration on all flow functions
pos_top = h(1).Position;
upper_bound = pos_top(2)+pos_top(4);
pos_bottom = h(3).Position;
lower_bound = pos_bottom(2);
combined_pos = [pos_top(1),lower_bound,pos_top(3),upper_bound-lower_bound];
combined_panel = axes('Units','Centimeters','Position',combined_pos);
plot(NaN,NaN,'r','LineWidth',3);
hold on
plot(NaN,NaN,'b','LineWidth',3);
plot(NaN,NaN,'Color',[0.5 0.5 0.5],'LineWidth',3);
xlim([-pi,pi]);
ylim([0 1]);
set(combined_panel,'XColor','None','YColor','None')
combined_panel.Color = 'None';
rectangle('Position',[-pi/2,0,pi,1],'FaceColor',[0 0 1 0.1],'EdgeColor','None')
bold_txts(1) = text(h(1), 0.125, 1.1, 'Retraction', 'Units', 'normalized', ...
     'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', 'FontSize', 12);
bold_txts(2) = text(h(1), 0.5, 1.1, 'Protraction', 'Units', 'normalized', ...
     'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', 'FontSize', 12);
bold_txts(3) = text(h(1), 0.875, 1.1, 'Retraction', 'Units', 'normalized', ...
     'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', 'FontSize', 12);
uistack(combined_panel,'top')
leg = legend({'Contact';'Non-contact';'Bootstrapped'},'Units','Centimeters','Color','None');
leg.NumColumns = 3;
set(leg,'Position',[5.3 26.1 3.255 0.5])



%% Cycle based one-dimensional (=alpha) visit rates:
side_pref_thresh = 0.9;
pre_cont_thresh_coeff = 1;

res = getDataPerCycle(allJobs,px2mm,side_pref_thresh,pre_cont_thresh_coeff,cont_thresh_px,extend_to_cycle_end_flag);

%% Cycle-based analysis

f_names = {'contact';'same_side';'other_side'};


clrs = cell(0);
clrs{1} = [1,0.65 0; 1, 0.35, 0; 1, 0, 0];
clrs{2} = [1, 0.8, 1; 1, 0.4, 1; 1, 0, 1];
clrs{3} = [0.8, 0.8, 1; 0.4, 0.4, 1; 0, 0, 1];

% PDFs on the units circle, collapsed across cycles, by cycle 'type' (first
% / interim / last)
fh = gobjects(3,3);

ttls = {'Contact';'Same side';'Other side'};
p_hist_bins = 73;
r_scale = 0.75;
all_hist_obs = cell(0);
polar_n_txt = [];
for k = 1:length(f_names)
    grp_cycles = res.(f_names{k});
    grp_cycles = grp_cycles(cellfun(@(x) length(x)>1,grp_cycles));

    first = cellfun(@(x) x{1},grp_cycles,'Uni',0);
    first_vals = cell2mat(cellfun(@(x) x(:)',first,'Uni',0));
    last = cellfun(@(x) x{end},grp_cycles,'Uni',0);
    last_vals = cell2mat(cellfun(@(x) x(:)',last,'Uni',0));
    interim = grp_cycles(cellfun(@(x) length(x)>2,grp_cycles));
    interim = cellfun(@(x) x(2:(end-1)),interim,'Uni',0);
    interim = [interim{:}];
    interim_vals = cell2mat(cellfun(@(x) x(:)',interim,'Uni',0));
    all_hist_obs(k,:) = {first_vals;interim_vals;last_vals}';

    fh(k,1) = axes('Units','centimeters');
    fh(k,1).Position = [1.5 13.5-3.5*(k-1) 2.75 2.75];
    polarhistogram(first_vals,-pi:(2*pi/p_hist_bins):pi,'Normalization','pdf','EdgeColor','None','FaceColor',clrs{k}(1,:));
    n_txt = ['n = ' regexprep(sprintf('%d', length(first_vals)), '\d(?=(\d{3})+$)', '$&,')];
    polar_n_txt(end+1) = text(gca,-pi/2,0.5,n_txt,'VerticalAlignment',...
        'middle','HorizontalAlignment','center','Color',[1 1 1]*0.5);
    if k == 1
        pos_to_neg_rad_func = @(x) x-(x > pi)*2*pi;
        theta_ticks = 0:pi/4:2*pi;
        theta_tick_labels = pos_to_neg_rad_func(theta_ticks);
        set(gca,'ThetaAxisUnits','radians','ThetaTick',theta_ticks,'ThetaTickLabel',makePolarTicks(theta_tick_labels),'TickLabelInterpreter','latex');
    else
        set(gca,'ThetaTickLabel','');
    end
    set(gca,'RLim',[0 r_scale]);
    if ~(k==1)
        set(gca,'RTickLabel','');
    end
    annot = annotation('textbox');
    set(annot,'Units','centimeters','Position',[1 13.5-3.75*(k-1)-0.5 4 0.5]);
    set(annot,'String',ttls{k},'EdgeColor', 'none', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle','Rotation',90,'FontSize',12);
    if k == 1
        bold_txts(end+1) = title('First','Units','centimeters','Position',[1.25 3.25 0]);
    end
    fh(k,2) = axes('Units','centimeters');
    fh(k,2).Position = [1.5+3.25 13.5-3.5*(k-1) 2.75  2.75];
    polarhistogram(interim_vals,-pi:(2*pi/p_hist_bins):pi,'Normalization','pdf','EdgeColor','None','FaceColor',clrs{k}(2,:));
    n_txt = ['n = ' regexprep(sprintf('%d', length(interim_vals)), '\d(?=(\d{3})+$)', '$&,')];
    polar_n_txt(end+1) = text(gca,-pi/2,0.5,n_txt,'VerticalAlignment',...
        'middle','HorizontalAlignment','center','Color',[1 1 1]*0.5);
    set(gca,'ThetaTickLabel','');
    if k == 1
        bold_txts(end+1) = title('Interim','Units','centimeters','Position',[1.25 3.25 0]);
    end
    set(gca,'RLim',[0 r_scale],'RTickLabel','');
    
    fh(k,3) = axes('Units','centimeters');
    fh(k,3).Position = [1.5+6.5 13.5-3.5*(k-1) 2.75 2.75];
    polarhistogram(last_vals,-pi:(2*pi/p_hist_bins):pi,'Normalization','pdf','EdgeColor','None','FaceColor',clrs{k}(3,:));
    n_txt = ['n = ' regexprep(sprintf('%d', length(last_vals)), '\d(?=(\d{3})+$)', '$&,')];
    polar_n_txt(end+1) = text(gca,-pi/2,0.5,n_txt,'VerticalAlignment',...
        'middle','HorizontalAlignment','center','Color',[1 1 1]*0.5-0.5*(k==1));
    set(gca,'ThetaTickLabel','');
    if k == 1
        bold_txts(end+1) = title('Last','Units','centimeters','Position',[1.25 3.25 0]);
    end
    set(gca,'RLim',[0 r_scale],'RTickLabel','');

end

% Repeat, but don't collapse - instead use on value per cycle and plot CDFs
cyc_func = @(x) dispersionMinimization(x(:));
cdf_ranges = 0:0.01:pi;


ttls = {'Contact';'Same-side';'Other-side'};
fhh = gobjects(1,3);
for k = 1:length(f_names)
    grp_cycles = res.(f_names{k});
    grp_cycles = grp_cycles(cellfun(@(x) length(x)>1,res.(f_names{1})));

    first = cellfun(@(x) x{1},grp_cycles,'Uni',0);
    interim = grp_cycles(cellfun(@(x) length(x)>2,grp_cycles));
    interim = cellfun(@(x) x(2:(end-1)),interim,'Uni',0);
    interim = [interim{:}];
    last = cellfun(@(x) x{end},grp_cycles,'Uni',0);

    first_to_last = [{first};{interim}; {last}];
    cdf_range = cdf_ranges;
    for mm = 1:length(first_to_last)
        cdf_vals = cellfun(@(x) cyc_func(x),first_to_last{mm});
        all_cycles{k,mm} = cdf_vals;
        if mm == 1
            fh(1,k) = axes('Units','centimeters');
        end
        fh(1,k).Position = [13 13.75-3.75*(k-1) 2.75 2.75];
        axes(fh(1,k));
        hold on
        [y,x] = histcounts(cdf_vals,cdf_range,'Normalization','cdf');
        x = x(1:end-1)+0.5*mode(diff(x));
        plot(x,y,'Color',clrs{k}(mm,:));
        ylabel('CDF')
        x_ticks = linspace(0,pi,5);
        x_tick_lbls = makePolarTicks(x_ticks);
        x_ticks(3) = pi/2-alpha_star_displacement;
        x_tick_lbls{3} = '$\alpha^{*}$';
        set(gca,'XTick',x_ticks,'XTickLabel',x_tick_lbls,'TickLabelInterpreter','latex');
        plot(pi/2-alpha_star_displacement*[1,1],[0 1],'k--');
        axis([min(cdf_range) max(cdf_range) 0 1]);
        if k == 3
            xlabel('$\arg\min \mathrm{MAD}(\alpha)$', 'Interpreter', 'latex');
        end
    
    end
    N_first = num2str(length(all_cycles{k,1}));
    N_interim = num2str(length(all_cycles{k,2}));
    N_last = num2str(length(all_cycles{k,3}));
    leg = legend({['First (n = ',N_first,')'];'';...
                  ['Interim (n = ',N_interim,')'];''; ...
                  ['Last (n = ',N_last,')']},'Location','eastoutside');
    set(leg,'Units','Centimeters');
    set(leg,'Position',get(leg,'Position')+[1.55 -0.5 0 0]);

end

% Caclulate some statistics for reporting (Wilcoxon and Levene):
final_cyc_grps = all_cycles(:,3);
final_cyc_grp_labels = cellfun(@(x,y) y*ones(size(x)),final_cyc_grps,num2cell(1:3)','Uni',0);
final_cyc_grps = [final_cyc_grps{:}]';
final_cyc_grp_labels = [final_cyc_grp_labels{:}]';
grp1 = final_cyc_grps(final_cyc_grp_labels == 1);
grp2 = final_cyc_grps(final_cyc_grp_labels == 2);
grp3 = final_cyc_grps(final_cyc_grp_labels == 3);
grp1_grp2 = [grp1; grp2];
grp1_grp3 = [grp1; grp3];
group_12_labels = [ones(size(grp1)); 2*ones(size(grp2))];
group_13_labels = [ones(size(grp1)); 2*ones(size(grp3))];
[p_dual_var_12, tbl] = vartestn(grp1_grp2, group_12_labels, ...
    'TestType', 'BrownForsythe', 'Display', 'off');
test_stat_12 = tbl.fstat;
[p_dual_var_13, tbl] = vartestn(grp1_grp3, group_13_labels, ...
    'TestType', 'BrownForsythe', 'Display', 'off');
test_stat_13 = tbl.fstat;
[p_contact_vs_alpha_star,~,rank_stats] = signrank(all_cycles{1,3}-(pi/2-alpha_star_displacement));
w_val = sum(rank_stats.signedrank(rank_stats.signedrank > 0)); % Compute W
p_contact_vs_alpha_star = [p_contact_vs_alpha_star,w_val,length(all_cycles{1,3})];

%% Permutation test: last cycle by shape and texture

% Parse to relevant trials
grp_cycles = res.contact;
shape = res.shape;
texture = res.texture;
shape = shape(cellfun(@(x) length(x)>1,grp_cycles));
texture = texture(cellfun(@(x) length(x)>1,grp_cycles));
grp_cycles = grp_cycles(cellfun(@(x) length(x)>1,grp_cycles));
last = cellfun(@(x) x{end}',grp_cycles,'Uni',0); % Get last cycle

obj_grps = {last(shape == 0 & texture == 0);...
           last(shape == 0 & texture == 1);...
            last(shape == 1 & texture == 0)};
obj_grps_trials = cellfun(@(y) cellfun(@(x) dispersionMinimization(x),y),obj_grps,'Uni',0);
obj_grps = cellfun(@(x) cell2mat(x),obj_grps,'Uni',0);
obj_grps = [[obj_grps{:}];obj_grps];
obj_grps = obj_grps([1,4,2,3]);
obj_grps_trials = [[obj_grps_trials{:}];obj_grps_trials];
obj_grps_trials = obj_grps_trials([1,4,2,3]);

% Plot polar-distributions
clrs{4} = [1, 0, 0;  0.2, 1, 0.2; 0, 0.7, 0; 0, 0.3, 0];


ttls = {{'All';'Objects'};{'Smooth';'Rectangular'};{'Smooth';'Circular'};{'Rough';'Circular'}};
for obj_num = 1:(length(obj_grps)-1)
    ph(obj_num) = polaraxes('Units','centimeters','Position',[1.5+3.25*(obj_num-1) 1.5 2.75 2.75]);
    polarhistogram(obj_grps{obj_num},-pi:(2*pi/p_hist_bins):pi,'Normalization','pdf','EdgeColor','None','FaceColor',clrs{4}(obj_num,:));
    set(gca,'ThetaTickLabel','','Rlim',[0 r_scale]);
    if obj_num~=1
        set(ph(obj_num),'RTickLabel',[]);
    end
    N_text_hist = ['n = ',num2str(length(obj_grps{obj_num}))];
    n_txt = ['n = ' regexprep(sprintf('%d', length(obj_grps{obj_num})), '\d(?=(\d{3})+$)', '$&,')];
    polar_n_txt(end+1) = text(gca,-pi/2,0.5,n_txt,'VerticalAlignment',...
        'middle','HorizontalAlignment','center','Color',[1 1 1]*0.5*((obj_num==2)|(obj_num==4)));
 
    title(ttls{obj_num});

end


% Plot median-distributions
ph(end+1) = axes('Units','centimeters','Position',[13 1.5 3 3]);
leg_ttls = {{'All'};{'Smooth Rect.'};{'Smooth Cir.'};{'Rough Circ.'}};
leg_ttls = cellfun(@(x,y) [x{1},char(10),' (n = ',num2str(length(y)),')'],leg_ttls,obj_grps_trials,'Uni',0);
for obj_num = 1:(length(obj_grps_trials)-1)
    
    [y,x] = histcounts(obj_grps_trials{obj_num},cdf_range,'Normalization','cdf');
    x = x(1:end-1)+0.5*mode(diff(x));
    obj_distribs(obj_num) = plot(x,y,'Color',clrs{4}(obj_num,:));
    uistack(obj_distribs(obj_num),'bottom')
    hold on
end
x_ticks = linspace(0,pi,5);
x_tick_lbls = makePolarTicks(x_ticks);
x_ticks(3) = pi/2-alpha_star_displacement;
x_tick_lbls{3} = '$\alpha^{*}$';
plot(pi/2-alpha_star_displacement*[1,1],[0 1],'k--');
set(gca,'XTick',x_ticks,'XTickLabel',x_tick_lbls,'TickLabelInterpreter','latex');
leg = legend(obj_distribs,leg_ttls(1:(end-1)));
set(leg,'Units','centimeters','Position',[15.2 2.2 3.5 1]);
ylabel('CDF')
xlabel('$\arg\min \mathrm{MAD}(\alpha)$', 'Interpreter', 'latex');

% Permutation tests
if get_permutation_p_val_flag


% Simple Kruskal Wallis:
obj_group_trials_labels = cellfun(@(x) ones(size(x)),obj_grps_trials,'Uni',0);
obj_group_trials_labels = cellfun(@(x,y) x.*y,obj_group_trials_labels,...
    num2cell(0:(length(obj_grps_trials)-1))','Uni',0);
obj_group_trials_labels = cell2mat(obj_group_trials_labels');
obj_group_trials_labels = obj_group_trials_labels...
    (obj_group_trials_labels>0);

% Permutation for total deviations (Kruskal-Wallis-like):
rng(121987,'twister');
all_vals = obj_grps_trials(2:end)';
all_vals = [all_vals{:}];
grand_median = median(all_vals);
group_meds = cellfun(@median, obj_grps_trials(2:end))';
T_obs = mean(abs(group_meds - grand_median));
group_sizes = cellfun(@length, obj_grps_trials(2:end))';
perm_T = [];
for bb = 1:1000
    shuffled = all_vals(randperm(length(all_vals)));
    
    % Split into permuted groups
    idx = [0, cumsum(group_sizes)];
    perm_groups = arrayfun(@(j) shuffled(idx(j)+1:idx(j+1)), 1:3, 'Uni', 0);
    
    % Compute test statistic
    perm_meds = cellfun(@median, perm_groups);
    perm_T(bb) = mean(abs(perm_meds - median(shuffled)));
end
perm_p_kruskal = mean(perm_T >= T_obs);

%  Objects differing by shape:
perm_p_pairs = [];
pair_ids = [2,3];
sim_res = [];
for mm = 1:1000
    grp1 = obj_grps_trials{pair_ids(1)};
    grp2 = obj_grps_trials{pair_ids(2)};
    agg_grp = [grp1,grp2];
    shuff_idx = randperm(length(agg_grp));
    N1 = length(grp1);
    grp1 = agg_grp(shuff_idx(1:N1));
    grp2 = agg_grp(shuff_idx(N1+1:end));
    sim_res(mm) = abs(median(grp1)-median(grp2));
end
grp1 = obj_grps_trials{pair_ids(1)};
grp2 = obj_grps_trials{pair_ids(2)};
gt = abs(median(grp1)-median(grp2));
perm_p_pairs(end+1) = mean(sim_res>gt);


%  Repeat for texture, correcting for small group size:
perm_p_pairs_size_corr = [];
pair_ids = [3,4];
sim_res = [];
for mm = 1:1000
    grp1 = obj_grps_trials{pair_ids(1)};
    grp2 = obj_grps_trials{pair_ids(2)};
    grp1 = randsample(grp1,length(grp2),0);
    agg_grp = [grp1,grp2];
    shuff_idx = randperm(length(agg_grp));
    N1 = length(grp1);
    grp1 = agg_grp(shuff_idx(1:N1));
    grp2 = agg_grp(shuff_idx(N1+1:end));
    sim_res(mm) = abs(median(grp1)-median(grp2));
end
grp1 = obj_grps_trials{pair_ids(1)};
grp2 = obj_grps_trials{pair_ids(2)};
gt = abs(median(grp1)-median(grp2));
perm_p_pairs_size_corr(end+1) = mean(sim_res>gt);

end


%% Cleaning up text formatting

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
for k = 1:length(bold_txts)
    set(bold_txts(k),'FontWeight','bold')
end
drawnow expose

%% Report statistics
disp('Fig 3 stats:')
% Datapoints:
disp('Fig 3a/1b, N_cont = ');
disp(length(flow_data{1}));
disp('Fig 3a/1b, N_non_cont = ');
disp(length(flow_data{2}));
disp(' ');

% Cycle variances:
cont_cycles_N = cellfun(@(x) length(x),all_cycles(1,:));
same_side_cycles_N = cellfun(@(x) length(x),all_cycles(2,:));
other_side_cycles_N = cellfun(@(x) length(x),all_cycles(3,:));

disp('Fig 3d, Brown-Forsythe test, last cycle - cont vs same side');
disp(['N_cont = ',num2str(cont_cycles_N(end))]);
disp(['N_same_side = ',num2str(same_side_cycles_N(end))]);
disp(['F = ',num2str(test_stat_12)]);
disp(['p = ',num2str(p_dual_var_12)]);
disp(' ');
disp('Fig 3d, Brown-Forsythe test, last cycle - cont vs other side');
disp(['N_cont = ',num2str(cont_cycles_N(end))]);
disp(['N_other_side = ',num2str(other_side_cycles_N(end))]);
disp(['F = ',num2str(test_stat_13)]);
disp(['p = ',num2str(p_dual_var_13)]);

disp(' ');
disp('Is last contact centered around alpha-star? p =');
disp('Fig 3d, Wilcoxon''s signed-rank test, last cycle cont vs (parameter) alpha-star');
disp(['N_cont = ',num2str(cont_cycles_N(end))]);
disp(['W = ',num2str(p_contact_vs_alpha_star(2))]);
disp(['p = ',num2str(p_contact_vs_alpha_star(1))]);

if get_permutation_p_val_flag
    
    disp('Permutation test for mean of deviations vs null (per trial):');
    disp('Tested with 1000 iterations:');
    disp(['p = ',num2str(perm_p_kruskal)]);

    disp('Permutation test for median, smooth rect vs smooth circ:');
    disp('Tested with 1000 iterations:');
    disp(['p = ',num2str(perm_p_pairs)]);
    disp('Permutation test for median, rough circ vs smooth_circ:');
    disp('Tested with 1000 iterations:');
    disp('Controlling for rough circ small group size:')
    disp(['p = ',num2str(perm_p_pairs_size_corr)]);

    
end




end