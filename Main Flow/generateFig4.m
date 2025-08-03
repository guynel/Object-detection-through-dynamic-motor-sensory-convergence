function [allJobs] = generateFig4(allJobs,px2mm,extend_to_cycle_end_flag,cont_thresh_px,alpha_star_val)

%% 
first_and_last_flag = 0; % Do we also want to plot interim cycles?

%% Get cycles:
side_pref_thresh = 0.9;
pre_cont_thresh_coeff = 5;
res = getDataPerCycle(allJobs,px2mm,side_pref_thresh,pre_cont_thresh_coeff,cont_thresh_px,extend_to_cycle_end_flag);
% Extract relevant variables:
alpha = res.contact;
r = res.r.contact;
RD = res.RD.contact;
vars = {alpha; r; RD};


%% Function choice for processing data points into parameters
top_prctile = 90;
median_func = @(x) median(x(:));
prctile_func = @(x) prctile(x(:),top_prctile);
r_func = @(x) log(prctile_func(x));
RD_func = @(x) (median_func(x));


alpha_star_displacement = pi/2-alpha_star_val;
dispersion_func = @(x) median(min(...
                            abs(circ_dist(x(:),pi/2-alpha_star_displacement)),...
                            abs(circ_dist(x(:),-pi/2-alpha_star_displacement))));
funcs = {dispersion_func; r_func; RD_func};

%% Variable ranges for plotting

var_rngs = [[0,pi-(pi/2-alpha_star_displacement)]; [0, 15]; [0, 4]];
var_rngs = [[0,pi-(pi/2-alpha_star_displacement)]; [-1.5,3.5]; [0, 4]];
x_lbls_with_units = {{'$\mathrm{MAD}(\alpha; \alpha^{*})$';'[radians]'},{'log. radius $r$ [magnitude]';'[$\|(\dot{\theta}, \dot{\kappa})_{scaled}\|$]'},'RD [cm]'};


%% Plot first and second rows of Fig. 4
% Plot first/interim/last CDF panels (calls a designated function)
main_fig = makeFullPagePDF();
step_size = [0.08, 0.5, 0.2];
step_size = [0.025 0.05 0.05];
p_first_last_handles = cell(0);
for var_num = 1:length(vars)
    CDF_panels(var_num) = axes();
    set(CDF_panels(var_num),'Units','Centimeters','Position',[2+4.3*(var_num+(var_num~=1)-1) 21 2.6 2]);
    x_rng = var_rngs(var_num,1):step_size(var_num):var_rngs(var_num,end);
    [p_vals_first_last(var_num,:),p_first_last_handles(end+1,:)] = plotFirstToLastCDFs(vars{var_num},funcs{var_num},CDF_panels(var_num),x_rng,x_lbls_with_units{var_num},first_and_last_flag);
    if var_num == 1
        set(gca,'XTick',linspace(0,pi/2,5),'XTickLabel',makePolarTicks(linspace(0,pi/2,5)),'TickLabelInterpreter','latex');
        n_first = num2str(sum(cellfun(@(x) length(x)~=1,vars{var_num})));
        n_last = n_first;
        n_single = num2str(sum(cellfun(@(x) length(x)==1,vars{var_num})));
        leg = legend({['First cycles (n = ',n_first,')'];...
            ['Last cycles (n = ',n_last,')'];...
            ['Single-contact \newline cycles  (n = ',n_single,')']},'Location','northoutside');
        set(leg,'Units','centimeters','Position',[3 23.75 1 1]);

    end
    x_lbl_handle = get(gca,'XLabel');
    x_lbl_handle.Units = 'centimeters';
    x_lbl_handle.Position = x_lbl_handle.Position-[0 0.5 0];
end
leg.ItemTokenSize = [15,1];


% Plot cycle depende panels (calls a designated function)
step_size = [0.01 0.01 0.01];
mult_fact = [20 50 50];
% For the function MAD(alpha,alpha-star), data is close to zero, so we
% don't want to use a kernel without accounting for that; however kernel
% use is purely for visualization anyway; linear regression is calculated
% on raw data
ks_opts{1} = {'Support',[0,pi],'BoundaryCorrection','log'}; % Be a tiny bit more careful with alpha, as it does near it's value's boundaries
ks_opts{2} = {};
ks_opts{3} = {};
% Dependence on cycle num.
for var_num = 1:length(vars)
    CDF_panels_3(var_num) = axes();
    set(CDF_panels_3(var_num),'Units','Centimeters','Position',[2+4.3*(var_num+(var_num~=1)-1) 15.25 2.6 2]);
    x_rng = var_rngs(var_num,1):step_size(var_num):var_rngs(var_num,end);
    [p_vals_slopes(var_num,:),lin_reg_handles(var_num,:)] = plotCycleDependence(vars{var_num},x_rng,funcs{var_num},mult_fact(var_num),ks_opts{var_num});
    
    if var_num ~= 1
        xlabel('N_{cycles}')
    else
        xlabel('N_{cycles}');
        set(gca,'YTick',linspace(0,pi/2,5),'YTickLabel',makePolarTicks(linspace(0,pi/2,5)),'TickLabelInterpreter','latex');
        hh = findobj(gca, 'Type', 'patch');
        leg = legend(hh([6,1]),{'First cycles';'Last cycles'},'Location','northoutside');
        set(leg,'Units','centimeters','Position',[2.75 17.5 1 1]);
    end
    x_tick_lbls = cellfun(@(x) num2str(x),num2cell(1:5),'Uni',0);
    x_tick_lbls{end} = ['$\geq',x_tick_lbls{end},'$'];
    set(gca,'XTick',1:5,'XTickLabel',x_tick_lbls,'TickLabelInterpreter','latex');
    ylabel(x_lbls_with_units{var_num},'Interpreter','latex');
    
end

leg.ItemTokenSize = [10,3];

%% Brown forsythe to show that r's variability increases

r_data = vars{2};
r_multi_cycle = r_data(cellfun(@(x) length(x)>1,r_data));
r_first = cellfun(@(x) x{1},r_multi_cycle,'Uni',0); % First
r_last = cellfun(@(x) x{end},r_multi_cycle,'Uni',0); % Last
r_first = cellfun(@(x) r_func(x),r_first);
r_last = cellfun(@(x) r_func(x),r_last);
r_labels = [ones(size(r_first)),2*ones(size(r_last))];
[p_r_brownForsythe, bf_tbl]  = vartestn([r_first,r_last]', r_labels',...
    'TestType', 'BrownForsythe', 'Display', 'off');
F_r_brownForsythe = bf_tbl.fstat;

% While we're at it, get last RD's median, IQR:
RD_data = vars{3};
RD_multi_cycle = RD_data(cellfun(@(x) length(x)>1,RD_data));
RD_last = cellfun(@(x) x{end},RD_multi_cycle,'Uni',0); % Last
RD_last = cellfun(@(x) RD_func(x),RD_last);
last_RD_range = prctile(RD_last,[25,50,75]);

%% We soon want to do the same for approach velocity, let's extract head-variables
all_d_hc_y = cell(0);
warning('off')
for k = 1:length(res.trial_num)
    trial_num = res.trial_num(k);
    tt = (res.tt.contact{k});
    d_hc_vec_y = cell(0);
    hc = allJobs{trial_num}.Tracks.headCenter;
    hc = abs(hc-allJobs{trial_num}.stimContour.stimCenter'); % Stimulus position is zero
    try
        lp_hc = lowpassFilt(hc,'Fst',100);
        if ~all(isnan(lp_hc))
            hc = lp_hc;
        end
    end
    for kk = 1:size(tt,1)
        t_vals = tt(kk,:);
        t_vals = t_vals(1):t_vals(end);
        cyc_hc = hc(:,t_vals);
        
        cyc_hc = 0.1*cyc_hc/px2mm; % Transform to cm
        d_hc_vec_y{kk} = 500*(multiStepDiff(cyc_hc(2,:),1)); % Account fot sampling rate
    end
    all_d_hc_y{k} =  d_hc_vec_y;
end
warning('on')
head_vars{1} = all_d_hc_y;

%% Function choice for processing head
head_func{1} = @(x) nanmean(x);

%% Repeat above plotting - for approach velocity

x_lbl = {'Approach velocity [cm/s]'};
x_rng = -30:0.1:10;
for var_num = 1
    CDF_panels = axes();
    set(CDF_panels,'Units','Centimeters','Position',[2+4.3 21 2.6 2]);
    [p_vals_first_last(end+1,:),p_first_last_handles(end+1,:)] = plotFirstToLastCDFs(head_vars{var_num},head_func{var_num},CDF_panels,x_rng,x_lbl,first_and_last_flag);
    
    x_lbl_handle = get(gca,'XLabel');
    x_lbl_handle.Units = 'centimeters';
    x_lbl_handle.Position = x_lbl_handle.Position-[0 0.5 0];
    drawnow expose
end


x_lbl = {'Approach';'velocity [cm/s]'};
x_rng = -30:0.1:10;
mult_fact = 30;
% Dependence on cycle num.
ks_opts{4} = {};
for var_num = 1
    CDF_panels = axes();
    set(CDF_panels,'Units','Centimeters','Position',[2+4.3 15.25 2.6 2]);
    [p_vals_slopes(end+1,:),lin_reg_handles(end+1,:)] = plotCycleDependence(head_vars{var_num},x_rng,head_func{var_num},mult_fact,ks_opts{4});
    x_tick_lbls = cellfun(@(x) num2str(x),num2cell(1:5),'Uni',0);
    x_tick_lbls{end} = ['$\geq',x_tick_lbls{end},'$'];
    set(gca,'XTick',1:5,'XTickLabel',x_tick_lbls,'TickLabelInterpreter','latex');
    xlabel('N_{cycles}')
end
%set(gca,'XTickLabel',[]);
ylabel(x_lbl);


%% Get correlation among whisker-variables, and with head (Panels e,f)
% Whisker related variables, function labels
all_vars = [vars];
all_funcs = [funcs];
var_names = {'$\alpha$';'log($r$)';'RD'};

% Go over all pairs of whisker-related variables
pairs_id = nchoosek(1:3,2);

% Pre-allocate for keeping r-values
first_panel_r_vals = [];
asterisk_handles_first_panel = gobjects(0);
for k = 1:size(pairs_id,1)
    % Get variables' indices:
    first_var_idx = pairs_id(k,1);
    sec_var_idx = pairs_id(k,2);

    % Make axis for the pair of variables
    first_panel_axes(k) = axes('Units','centimeters','Position',[2.75+1.75*(k-1) 9.5 1.75 3]);
    pair_lbls = cell(0);

    for mm = 1:3 % Iterate over correlation, time-aligned, lagged, leading

        first_var = all_vars{first_var_idx}; % First variable
        sec_var = all_vars{sec_var_idx}; % Second variable

        if mm == 1 % Aligned
            first_vals = [first_var{:}]; 
            sec_vals = [sec_var{:}];

        elseif mm == 2 % Current first predicts future second (first leads)
            first_vals = first_var(cellfun(@(x) length(x)~=1,first_var));
            first_vals = cellfun(@(x) x(1:end-1),first_vals,'Uni',0);
            first_vals = [first_vals{:}];

            sec_vals = sec_var(cellfun(@(x) length(x)~=1,sec_var));
            sec_vals = cellfun(@(x) x(2:end),sec_vals,'Uni',0);
            sec_vals = [sec_vals{:}];

        elseif mm == 3 % Current second predicts future first (second leads)
            first_vals = first_var(cellfun(@(x) length(x)~=1,first_var));
            first_vals = cellfun(@(x) x(2:end),first_vals,'Uni',0);
            first_vals = [first_vals{:}];

            sec_vals = sec_var(cellfun(@(x) length(x)~=1,sec_var));
            sec_vals = cellfun(@(x) x(1:end-1),sec_vals,'Uni',0);
            sec_vals = [sec_vals{:}];

        end

        % Extract values
        first_vals = cellfun(@(x) all_funcs{first_var_idx}(x),first_vals);
        sec_vals = cellfun(@(x) all_funcs{sec_var_idx}(x),sec_vals);

        % Calculate correlation
        [cc,p_val_mat] = corrcoef(first_vals,sec_vals);
        r_val = abs(cc(2));
        first_panel_r_vals(end+1,:) = [r_val,p_val_mat(2),length(first_vals)];
        
        % Add relevant labels
        switch mm
            case 1
                pair_lbl = 'aligned';
            case 2
                pair_lbl = [var_names{first_var_idx}, ' leads'];
            otherwise
                pair_lbl = [var_names{sec_var_idx}, ' leads'];
        end
        % Plot
        pair_lbls{end+1} = pair_lbl;
        bar(mm,r_val,'FaceColor','k','EdgeColor','k','FaceAlpha',0.5+0.5*(mm~=1))
        hold on
        % Place asterisk (we'll remove if not significant)
        asterisk_handles_first_panel(end+1) = text(mm,r_val-0.05, '*', 'Color','w','FontSize', 16, 'HorizontalAlignment', 'center');
    end
    
    % Tidy:
    set(gca,'XTick',1:3,'XTickLabel',pair_lbls,'TickLabelInterpreter','Latex');
    ylim([0 0.7])
    if k ~= 1
        set(gca,'YColor','none','box','off')
    else
        set(gca,'box','off','XColor','k','YColor','k')
        ylabel('Correlation')
    end
    title(['(',var_names{first_var_idx},', ',var_names{sec_var_idx},')'],'Interpreter','latex');
end

% Re-arrange order, next plot becomes nicer:
all_vars = all_vars([1,3,2]);
all_funcs = all_funcs([1,3,2]);
var_names = var_names([1,3,2]);

% Now each whisker-related variable with approach-velocity
head_vars_idx = 1;
tick_lbls = {'aligned';'head lags';'head leads'};
sec_panel_r_vals = [];
axes('Units','centimeters','Position',[11.5 9.5 5.25 3]);
asterisk_handles_sec_panel = gobjects(0);
for k = 1:length(vars) % Iterate over whisker-relate variables
    first_var_idx = k; % Whisker-related
    sec_var_idx = head_vars_idx; % Approach-velocity
    
    for mm = 1:3
        % Same as before:
        first_var = all_vars{first_var_idx};
        sec_var = head_vars{sec_var_idx};

        if mm == 1 % Aligned
            first_vals = [first_var{:}];
            sec_vals = [sec_var{:}];

        elseif mm == 2 % Current whisker predict future head
            first_vals = first_var(cellfun(@(x) length(x)~=1,first_var));
            first_vals = cellfun(@(x) x(1:end-1),first_vals,'Uni',0);
            first_vals = [first_vals{:}];

            sec_vals = sec_var(cellfun(@(x) length(x)~=1,sec_var));
            sec_vals = cellfun(@(x) x(2:end),sec_vals,'Uni',0);
            sec_vals = [sec_vals{:}];

        elseif mm == 3 % Current head predict future whisker
            first_vals = first_var(cellfun(@(x) length(x)~=1,first_var));
            first_vals = cellfun(@(x) x(2:end),first_vals,'Uni',0);
            first_vals = [first_vals{:}];

            sec_vals = sec_var(cellfun(@(x) length(x)~=1,sec_var));
            sec_vals = cellfun(@(x) x(1:end-1),sec_vals,'Uni',0);
            sec_vals = [sec_vals{:}];

        end
        % Extract values:
        first_vals = cellfun(@(x) all_funcs{first_var_idx}(x),first_vals);
        sec_vals = cellfun(@(x) head_func{sec_var_idx}(x),sec_vals);
        % Calculate correlation:
        [cc,p_val_mat] = corrcoef(first_vals,sec_vals);
        r_val = abs(cc(2));
        
        % Get relevant labels:
        sec_panel_r_vals(end+1,:) = [r_val,p_val_mat(2),length(first_vals)];
        switch mm
            case 1
                pair_lbl = 'aligned';
            case 2
                pair_lbl = [var_names{first_var_idx}, ' leads'];
            otherwise
                pair_lbl = [var_names{sec_var_idx}, ' leads'];
        end
        % Plot
        bar(mm+(k-1)*4,r_val,'FaceColor','k','EdgeColor','k','FaceAlpha',0.5+0.5*(mm~=2));
        hold on
        asterisk_handles_sec_panel(end+1) = text(mm+(k-1)*4,r_val-0.05, '*', 'Color','w','FontSize', 16, 'HorizontalAlignment', 'center');
    end
    
    ylim([0 0.7])
    
    
   
end
% Previous panel (4e) was three-panel disguised as one; now its a single
% panel, so need to work harder to get x-ticks correct:
x_tick = ((1:3)+((1:3)'-1)*4)';
x_tick = x_tick(:);
x_tick_lbls = repmat(tick_lbls,1,3);
x_tick_lbls = [x_tick_lbls(:)]';
set(gca,'XTick',x_tick,'XTickLabel',x_tick_lbls,'TickLabelInterpreter','Latex');
ylabel('Correlation')
txt_pos_x = [1,2.6,4.4];
for k = 1:3
    txt = text('Units','centimeters','HorizontalAlignment','center','Position',[txt_pos_x(k) 3.5]);
    txt.String = var_names{k};
    txt.Interpreter = 'latex';
end
box off

% As the middle correlation (whisker-leads / head-lagged) peaks for all
% whisker related variables, check to see if they differ significantly
% Uses Fisher's z-test
p_fisher_head = [];
[~,p_fisher_head(end+1)] = fisher_z_test(sec_panel_r_vals(2,1),sec_panel_r_vals(5,1),sec_panel_r_vals(2,3),sec_panel_r_vals(5,3));
[~,p_fisher_head(end+1)] = fisher_z_test(sec_panel_r_vals(2,1),sec_panel_r_vals(8,1),sec_panel_r_vals(2,3),sec_panel_r_vals(8,3));
orig_p_fisher_head = p_fisher_head;
p_fisher_head = fdr_bh(p_fisher_head);
black_bar_asterisks = gobjects(0);
for k = 1:length(p_fisher_head)
    if p_fisher_head(k)
        plot([2,2+k*4],0.025+0.075*(k-1)+[1,1]*max(sec_panel_r_vals([2,5,8])),'k','LineWidth',2);
        hold on
        black_bar_asterisks(end+1) = text(mean([2,2+k*4]),0.045+0.085*(k-1)+max(sec_panel_r_vals([2,5,8])), '*', 'Color','k','FontSize', 16, 'HorizontalAlignment', 'center');
    end
end

%% Tidy up:
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

%% Significance tidying:
% Make asterisks correct size again
for k = 1:length(asterisk_handles_first_panel)
    asterisk_handles_first_panel(k).FontSize = 16;
end
for k = 1:length(asterisk_handles_sec_panel)
    asterisk_handles_sec_panel(k).FontSize = 16;
end
for k = 1:numel(p_first_last_handles)
    set(p_first_last_handles{k}(3),'FontSize',14,'FontWeight','Normal');
end
for k = 1:length(black_bar_asterisks)
    black_bar_asterisks(k).FontSize = 16;
end

% Delete asterisks from panel e if they don't survive FDR correction
delete(asterisk_handles_first_panel(~fdr_bh(first_panel_r_vals(:,2))'))
% Delete asterisks from panel f if they don't survive FDR correction
delete(asterisk_handles_sec_panel(~fdr_bh(sec_panel_r_vals(:,2))'))

%% Sorting p-vals
% The plotting order has approach-velocity at the end, but it should be
% second:
p_vals_first_last = p_vals_first_last([1,4,2,3],:);
p_first_last_handles = p_first_last_handles([1,4,2,3],:);
p_vals_slopes = p_vals_slopes([1,4,2,3],:);
fdr_p_vals_first_last = fdr_bh(p_vals_first_last);
fdr_p_vals_slopes = fdr_bh(p_vals_slopes);
lin_reg_handles = lin_reg_handles([1,4,2,3],:);

% Remove linear-regressions that don't survive FDR
delete(lin_reg_handles(~fdr_p_vals_slopes));

% Remove Wilcoxon (top row panels a-d) that don't survive FDR
for k = 1:numel(fdr_p_vals_first_last)
    if ~fdr_p_vals_first_last(k)
        for kk = 1:length(p_first_last_handles{k})
            delete(p_first_last_handles{k}(kk))
        end
    end
end

%% Report statistics:
disp('Fig. 4 stats:')

disp(' ');
helper_ttls = {'MAD(alpha)';'Approach Velocity';'radius';'RD'};
disp('Fig. 4a-b: single vs first and vs last')
for k = 1:4
    disp([helper_ttls{k},': ']);
    disp(['Single vs first: p = ',num2str(p_vals_first_last(k,1)),'; FDR_sig = ',num2str(fdr_p_vals_first_last(k,1))]);
    disp(['Single vs last: p = ',num2str(p_vals_first_last(k,2)),'; FDR_sig = ',num2str(fdr_p_vals_first_last(k,2))])
end


disp(' ');
helper_ttls = {'MAD(alpha)';'Approach Velocity';'radius';'RD'};
disp('Fig. 4a-b: slopes')
for k = 1:4
    disp([helper_ttls{k},': ']);
    disp(['Slope first: p = ',num2str(p_vals_slopes(k,1)),'; FDR_sig = ',num2str(fdr_p_vals_first_last(k,1))]);
    disp(['Slope last : p = ',num2str(p_vals_slopes(k,2)),'; FDR_sig = ',num2str(fdr_p_vals_first_last(k,2))])
end

disp('Fig. 4c, first vs last variance (Brown Forsythe):');
disp([' p = ',num2str(p_r_brownForsythe),';' ...
    ' F = ',num2str(F_r_brownForsythe)]);

disp('Fig. 4d, last cycles RD:')
disp(['median = ',num2str(last_RD_range(2)),'; IQR = ',...
    num2str(last_RD_range(1)),'-',num2str(last_RD_range(end))]);

disp(' ');
disp('Fig. 4e: r_values by order')
for k = 1:size(first_panel_r_vals,1)
    disp(['r = ',num2str(first_panel_r_vals(k,1)),'; p = ',num2str(first_panel_r_vals(k,2))]);
end
disp('Asterisks denote correlation significant after FDR correction for all 9 bars');

disp(' ');
disp('Fig. 4f: r_values by order')
for k = 1:size(sec_panel_r_vals,1)
    disp(['r = ',num2str(sec_panel_r_vals(k,1)),'; p = ',num2str(sec_panel_r_vals(k,2))]);
end
disp('Asterisks denote correlation significant after FDR correction for all 9 bars');

disp('Fig. 4f: head lagged across vairables ');
helper_ttls = {'alpha vs RD';'alpha vs r'};
for k = 1:length(orig_p_fisher_head)
    disp([helper_ttls{k},'; p = ',num2str(orig_p_fisher_head(k)),'; FDR_sig = ',num2str(p_fisher_head(k))])
end
end