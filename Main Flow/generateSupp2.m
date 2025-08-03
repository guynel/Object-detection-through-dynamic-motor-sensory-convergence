function [] = generateSupp2(allJobs,px2mm,extend_to_cycle_end_flag,cont_thresh_px,alpha_star_val)
% Extended Fig. 2

% Get cycles:
side_pref_thresh = 0.9;
pre_cont_thresh_coeff = 5;
res = getDataPerCycle(allJobs,px2mm,side_pref_thresh,pre_cont_thresh_coeff,cont_thresh_px,extend_to_cycle_end_flag);

%%
% Some parameters
NN = 200; % Moving window size
max_step = pi./2; % Same as in fig. 2
step_N = 3; % Use larger delta than fig. 2, for additional de-noising


trial_num_vec = res.trial_num;
side_num_vec = res.side_num;

% Collect data
all_seq = cell(0);
all_non_seq = cell(0);
for trial_num_idx = 1:length(trial_num_vec)
    trial_num = trial_num_vec(trial_num_idx);
    side_num = side_num_vec(trial_num_idx);
    trim_idx = allJobs{trial_num}.Tracks.trim_idx;
    if isnan(trim_idx)
        continue
    end
    dd = allJobs{trial_num}.Tracks.whiskerObjDist(:,1:trim_idx)<cont_thresh_px;
    side = allJobs{trial_num}.Tracks.side;
    side = mode(side,2);
    cont = dd./dd;
    cont_side = side == side_num;
    cont_side = cont_side./cont_side;
    cont = cont.*side;
    non_cont = (~dd)./(~dd);

    d_ang = allJobs{trial_num}.Tracks.diffAngles(:,1:trim_idx);
    d_curv = allJobs{trial_num}.Tracks.diffCurv(:,1:trim_idx);

    [alpha_vals,r_vals] = cart2pol(d_ang,d_curv);
    cont_alpha = alpha_vals.*cont;
    non_cont_alpha = alpha_vals.*non_cont;

    cont_side_t = res.tt.contact{trial_num_idx};
    cont_side_t(end) = trim_idx;
    
    % Timeseries pairs for cont:
    cont_cyc_seq = cell(0);
    for cyc_num = size(cont_side_t,1):-1:1
        tt = cont_side_t(cyc_num,1):cont_side_t(cyc_num,2);
        cont_alpha_cyc = cont_alpha(:,tt);
        d_alpha = cont_alpha_cyc(:,step_N:end);
        origin_alpha = cont_alpha_cyc(:,1:(end-(step_N-1)));
        d_alpha = [origin_alpha(:),d_alpha(:)];
        d_alpha = d_alpha(~any(isnan(d_alpha),2),:);
        cont_cyc_seq{end+1} = d_alpha;
    end
    all_seq{end+1} = cont_cyc_seq;

    % Timeseries pairs for non-cont:
    non_cont_cyc_seq = cell(0);
    for cyc_num = size(cont_side_t,1):-1:1
        tt = cont_side_t(cyc_num,1):cont_side_t(cyc_num,2);
        non_cont_alpha_cyc = non_cont_alpha(:,tt);
        d_alpha = non_cont_alpha_cyc(:,step_N:end);
        origin_alpha = non_cont_alpha_cyc(:,1:(end-(step_N-1)));
        d_alpha = [origin_alpha(:),d_alpha(:)];
        d_alpha = d_alpha(~any(isnan(d_alpha),2),:);
        non_cont_cyc_seq{end+1} = d_alpha;
    end
    all_non_seq{end+1} = non_cont_cyc_seq;
end

% Remove trials with single contact
rmv_idx = cellfun(@(x) length(x) == 1,all_seq);
all_seq(rmv_idx) = [];
all_non_seq(rmv_idx) = [];

% Aggregate everything, for cont; reversed so last contact is first; no
% more than 3 cycles
for kk = 1:length(all_seq)
    trial_seq = all_seq{kk};
    first_idx = find(cellfun(@(x) ~isempty(x),trial_seq),1);
    last_idx = find(cellfun(@(x) ~isempty(x),trial_seq),1,'last');
    trial_seq = trial_seq(first_idx:last_idx);
    trial_seq = trial_seq(1:min(length(trial_seq),3));
    trial_seq = [trial_seq,repmat({NaN(1,2)},1,3-length(trial_seq))];
    all_seq{kk} = trial_seq;
end


% Same for non-contact
for kk = 1:length(all_seq)
    trial_seq = all_non_seq{kk};
    first_idx = find(cellfun(@(x) ~isempty(x),trial_seq),1);
    last_idx = find(cellfun(@(x) ~isempty(x),trial_seq),1,'last');
    trial_seq = trial_seq(first_idx:last_idx);
    trial_seq = trial_seq(1:min(length(trial_seq),3));
    trial_seq = [trial_seq,repmat({NaN(1,2)},1,3-length(trial_seq))];
    all_non_seq{kk} = trial_seq;
end

% Make into vector, by cycle number
cont_cyc = cell(0);
non_cont_cyc = cell(0);
for k = 1:3 % Iterate over cycles
    temp = cellfun(@(x) x{k}', all_seq, 'UniformOutput', false);
    temp = [temp{:}]';
    temp(any(isnan(temp),2),:) = []; % Remove NaN
    cont_cyc{k} = temp;

    % Same for non-contact
    temp = cellfun(@(x) x{k}', all_non_seq, 'UniformOutput', false);
    temp = [temp{:}]';
    temp(any(isnan(temp),2),:) = [];
    non_cont_cyc{k} = temp;
end

% Plotting
fig = makeFullPagePDF();
ttls = {'Cycle = Last';'Cycle = -1';'Cycle = -2'};
rng(121987,'twister');
for cyc_num = 1:3
    % Get cycle data
    all_cont_d_alpha = cont_cyc{cyc_num};
    all_non_cont_d_alpha = non_cont_cyc{cyc_num};

    % Get everything set up:
    % Slightly different flow function here:
    flow_func = @(x,y) movmean(abs(circ_dist(x(:,2),y))<abs(circ_dist(x(:,1),y)),NN);
    % All data points:
    flow_data = {all_cont_d_alpha;all_non_cont_d_alpha};
    data_vec = [flow_data{1};flow_data{2}];
    % Corresponding labels (1 for contact):
    labels = [ones(size(flow_data{1},1),1);zeros(size(flow_data{2},1),1)];
    % Remove large steps:
    labels(abs(circ_dist(data_vec(:,1),data_vec(:,2)))>max_step) = [];
    data_vec(abs(circ_dist(data_vec(:,1),data_vec(:,2)))>max_step,:) = [];
    % How many contact pts per cycle
    output_Ns(cyc_num) = sum(labels == 1);

    % Plotting
    alpha_bounds = [0,pi]; % We only look at protraction
    axes('Units','Centimeters','Position',[15.5-5*(cyc_num-1) 10 3.5 2]);
    for dd = 1:(size(alpha_bounds,2)-1)
        % Remove data out of bounds"
        bin_idx = (data_vec(:,1)>alpha_bounds(dd)) & (data_vec(:,1)<alpha_bounds(dd+1));
        bin_vec = data_vec(bin_idx,:);
        bin_labels = labels(bin_idx);

        % Apply function to actual data:
        gt = bin_vec(logical(bin_labels),:); % Actual data, source + step
        [~,idx] = sort(gt(:,1)); % Sort (source and destination) by source
        gt = gt(idx,:);
        x_val = movmean(gt(:,1),NN); % Smooth the sources
        y_val = flow_func(gt,alpha_star_val); % Get flows
        plot(x_val,y_val,'r.') % Plot

        % Repeat for label-shuffled data
        hold on
        sim_res = [];
        for iter = 1:100
            % Out of entire population (cont+no_cont), randomly sample same
            % sized group:
            idx = randsample(length(bin_vec),sum(bin_labels),1);
            gt = bin_vec(idx,:);
            % Repeat procedure from above to get flow:
            [~,idx] = sort(gt(:,1));
            gt = gt(idx,:);
            x_val = movmean(gt(:,1),NN);
            y_val = flow_func(gt,alpha_star_val);
            q = plot(x_val,y_val,'Color',[0 0 0 0.05]);
            uistack(q,'down');
        end
    end

    % Sort out plot
    ylim([0 1])
    plot(pi/2*sign(alpha_star_val).*[1,1],ylim,'k--')
    plot(alpha_star_val.*[1,1],ylim,'r--')
    plot(xlim,[1,1]*0.5,'k--')

    % Sort out ticks, labels, etc
    drawnow
    x_ticks = linspace(0,pi,5);
    x_tick_lbls = makePolarTicks(x_ticks);
    x_ticks = [x_ticks(1:2),alpha_star_val,x_ticks(4:end)];
    x_tick_lbls = [x_tick_lbls(1:2),'$\alpha^*$',x_tick_lbls(4:end)];
    set(gca,'XTick',x_ticks,'XTickLabels',x_tick_lbls,'TickLabelInterpreter','latex')
    if cyc_num == 3
        ylabel('p(moving towards $\alpha^{*})$','Interpreter','latex')
        
    end
    xlabel('$\alpha$','Interpreter','latex');
    set(gca,'TickLabelInterpreter','latex')
    title(ttls{cyc_num});
    ylim([0.2 0.7])
end

% Tidy up:
child_list = get(gcf,'children');
for k = 1:length(child_list)
    try
        set(child_list(k),'FontSize',9,'FontWeight','Normal');
    end
end
txt_list = findall(gcf,'Type','Text');
for k = 1:length(txt_list)
    set(txt_list(k),'FontSize',9,'FontWeight','Normal');
end
axes_list = findall(gcf,'Type','Axes');
for k = 1:length(axes_list)
    set(axes_list(k),'FontSize',9,'FontWeight','Normal');
    set(get(axes_list(k),'XLabel'),'FontSize',9,'FontWeight','Normal');
    set(get(axes_list(k),'YLabel'),'FontSize',9,'FontWeight','Normal');
end
clrbar_list = findall(gcf,'Type','Colorbar');
for k = 1:length(clrbar_list)
    set(clrbar_list(k),'FontSize',9,'FontWeight','Normal');
end
annot_list = findall(gcf,'Type','textboxshape');
for k = 1:length(annot_list)
    set(annot_list(k),'FontSize',9,'FontWeight','Normal');
end
drawnow expose

disp('')
disp('Extended Fig. 2 stats:')
disp('Sample sizes, from last in reverse:');
disp(num2str(output_Ns))
end