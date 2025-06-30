function [allJobs] = generateFig5(allJobs,px2mm,pump_struct,extend_to_cycle_end_flag,alpha_star_displacement,cont_thresh_px)

%% Get cycles by type
full_cycle_idx = [pump_struct.complete_cycle];
double_pump_idx = [pump_struct.is_double_pump];
contact_idx = [pump_struct.is_cont];
tip_idx = [pump_struct.is_tip];
during_approach_idx = [pump_struct.during_approach];

single_free_air_idx = ~double_pump_idx & ~contact_idx & full_cycle_idx & during_approach_idx;
single_contact_idx = ~double_pump_idx & contact_idx & full_cycle_idx & during_approach_idx;
double_free_air_idx = double_pump_idx & ~contact_idx & ~tip_idx & full_cycle_idx & during_approach_idx;
touch_induced_pump_idx = double_pump_idx & contact_idx & tip_idx & full_cycle_idx & during_approach_idx;

%% Plot select examples
% In order: free air, contact, TIPs, TOPs
select_cycles_idx = [13,76,85,4,2277,1471, 1792,2132,2336,621,2202, 2099];
counter = 0;

fig = makeFullPagePDF();
max_t = 0;
ttls = {'Free-air';'Double-pump';'Contact';'TIP'};
for pump_idx = 1:length(select_cycles_idx)
    pump_num = select_cycles_idx(pump_idx);

    counter = counter+1;
    if rem(counter,3) == 1
        fh(ceil((1:12)/3)) = axes('Units','centimeters','Position',[1+0.8*floor(counter/3) 22 0.8 2.5]);
        ttl = title(ttls{ceil(counter/3)},'Rotation',45);
        ttl.Units = 'centimeters';
        ttl.HorizontalAlignment = 'left';
    end
    hold on

    t = pump_struct(pump_num).start_end;
    t = t(1):t(end);
    ang = abs(pump_struct(pump_num).lp_med_ang(t));
    ang = zscore(ang); % Get them in the same range, we only care about shape here
    [max_val,max_idx] = max(ang);
    ang = ang-max_val+3.25*(1+rem(counter,3));
    tt = (t-t(1));
    plot(tt,ang,'k');
    hold on
    dd = min(pump_struct(pump_num).wh_obj_dist(:,t));
    cont = dd<12;
    set(fh(ceil((1:12)/3)),'XTickLabel',[],'YTickLabel',[]);
    if sum(cont) == 1
        plot(tt,ang.*(cont./cont),'o-','LineWidth',2);
    else
        plot(tt,ang.*(cont./cont),'r-','LineWidth',2);
    end
    drawnow
    if counter == 1
        x_lbl = xlabel('Time');
        set(x_lbl,'Units','centimeters');
        set(x_lbl,'Position',get(x_lbl,'Position')+[1.5,0,0])
        ylabel('\theta (normalized)')
    end

    if ~rem(counter,3)
        xlim([-10 max(xlim).*1.1])
    end
end

%% Characteristic duration and amplitude
% For all bilateral cycles during approach, get duration and amplitude:
conds = [single_free_air_idx; single_contact_idx; double_free_air_idx; touch_induced_pump_idx];
ln_types = {'k-';'r-';'k-o';'r-o'};
all_cyc_amplitude = cell(0);
all_cyc_duration = cell(0);
all_cyc_duration_double = cell(0);
for cond_num = 1:size(conds,1)
    cond = conds(cond_num,:);


    cyc_amplitude = arrayfun(@(x) range(x.lp_med_ang), pump_struct(cond));
    if cond_num == 1
        fh(5) = axes('Units','Centimeters','Position',[6 22 2.5 2.5]);
        hold on
        plot(NaN,NaN,'r-');
        plot(NaN,NaN,'k-');
        plot(NaN,NaN,'ko');
    else
        axes(fh(5));
    end
    all_cyc_amplitude{end+1} = cyc_amplitude;
    [y,x] = histcounts(cyc_amplitude,0:2:100,'Normalization','cdf');
    x = x(1:end-1)+(0.5*mode(diff(x)));
    plot(x,y,ln_types{cond_num},'MarkerIndices',1:4:length(y));
    hold on
    xlabel('Cycle Amplitude [Deg.]');
    ylabel('CDF');

    cyc_duration = 2*arrayfun(@(x) diff(x.start_end), pump_struct(cond));
    if cond_num == 1
        fh(6) = axes('Units','Centimeters','Position',[10.5 22 2.5 2.5]);
        hold on
        plot(NaN,NaN,'k--');
    else
        axes(fh(6));
    end
    hold on
    all_cyc_duration{end+1} = cyc_duration;
    [y,x] = histcounts(cyc_duration,0:3:300,'Normalization','cdf');
    x = x(1:end-1)+(0.5*mode(diff(x)));
    plot(x,y,ln_types{cond_num},'MarkerIndices',1:8:length(y));
    hold on

    if cond_num == 1 | cond_num == 2
        cyc_duration = cyc_duration+cyc_duration';
        cyc_duration = cyc_duration(:)';
        all_cyc_duration_double{end+1} = cyc_duration;
        axes(fh(6));
        [y,x] = histcounts(cyc_duration,0:3:300,'Normalization','cdf');
        x = x(1:end-1)+(0.5*mode(diff(x)));
        plot(x,y,[ln_types{cond_num}(1),'--']);
        hold on
    end
    xlabel('Cycle Duration [ms]');
    ylabel('CDF');
end

% First col is free air, second col cont; first row single, second col
% double-pump:
all_cyc_amplitude = reshape(all_cyc_amplitude,2,2);
all_cyc_duration = reshape(all_cyc_duration,2,2);
% Is TIP longer-duration than contact?
[p_val_tip_vs_cont_dur,~,stats_tip_vs_cont_dur] = ranksum(all_cyc_duration{2,1},all_cyc_duration{2,2},'tail','left');
% By how much, on average?
TIP_elongated_prct = mean(all_cyc_duration{2,2})./mean(all_cyc_duration{2,1});

% Add legend
axes(fh(6))
leg = legend('Location','northoutside');
elements =  get(gca, 'Children');
leg_element_lbls = {'Free-air';'Contact';'Double-pump \newline free-air';'TIP';'Two free-air \newline cycles';'Two contact \newline cycles'};
leg = legend(elements([6,4,2,1,7,3]),leg_element_lbls,'Location','northoutside');
set(leg,'Units','centimeters','Position',[5.9 25.75 8 1])
leg.NumColumns = 3;

%% Work on the unilateral contact cycles, same as in previous figure
side_pref_thresh = 0.9;
pre_cont_thresh_coeff = 5;
res = getDataPerCycle(allJobs,px2mm,side_pref_thresh,pre_cont_thresh_coeff,cont_thresh_px,extend_to_cycle_end_flag);

%% When do TIPs occur temporally?

trial_num = [res.trial_num];
side_num = [res.side_num];
tt = [res.tt.contact];

trial_num_vec = [pump_struct.trial_num];
side_num_vec = [pump_struct.side_num];
t_vec = [pump_struct.start_end];
t_vec = reshape(t_vec,2,length(t_vec)/2)';

% Order all contact cycles by last, backwards
is_tip_vec = cell(0);
for k = 1:length(trial_num)
    is_tip_temp = [];
    for kk = 1:size(tt{k},1)
        trial_cond = trial_num_vec == trial_num(k);
        side_cond = side_num_vec == side_num(k);
        t_cond = all(t_vec == tt{k}(kk,:),2);
        idx = find(trial_cond & side_cond & t_cond');
        is_tip_temp(end+1) = pump_struct(idx).is_tip;
    end
    is_tip_vec{k} = is_tip_temp;
end
unilaterall_approach_TIPs = sum(cellfun(@(x) sum(x),is_tip_vec));

tip_ordinal = cellfun(@(x) -1*(length(x)-(1:length(x))),is_tip_vec,'Uni',0);
tip_ordinal_vec = [tip_ordinal{:}];
is_tip_vec = [is_tip_vec{:}];

% Extract and plot
fh(7) = axes('Units','Centimeters','Position',[14.75 22 2.5 2.5]);
bar_xy = [];
tip_bars_Ns = zeros(1,5);
for k = 0:-1:-4
    if k ~= -4
        p_tip = mean(is_tip_vec(tip_ordinal_vec == k));
        tip_bars_Ns(abs(k)+1) = sum(tip_ordinal_vec == k);
    else
        p_tip = mean(is_tip_vec(tip_ordinal_vec <= k));
        tip_bars_Ns(abs(k)+1) = sum(tip_ordinal_vec <= k);
    end
    bar(k,p_tip,'FaceColor',[1 1 1]*0.75);
    hold on
    bar_xy(end+1,:) = [k,p_tip];
end
set(gca,'XTick',[-4:1:0],'XTickLabel',{'\leq{-4}';'-3';'-2';'-1';'Last'})
xlabel('Cycle Num.')
ylabel('p(TIP)')


% Is last different from non-last? (using binomial)
grp1 = is_tip_vec(tip_ordinal_vec == 0);
grp2 = is_tip_vec(tip_ordinal_vec ~= 0);
table_res = [sum(grp1),sum(grp2);sum(~grp1),sum(~grp2)];
[~,p_binomial_last_cycle] = fishertest(table_res,'tail','right');
p_binomial_last_cycle = [p_binomial_last_cycle,mean(grp1),length(grp1),mean(grp2),length(grp2)];

% Is a logistic model better than a null (constant)?
tbl = table(bar_xy(:,1), bar_xy(:,2), 'VariableNames', {'x', 'y'});
mdl_full = fitglm(tbl, 'y ~ x', 'Distribution', 'binomial');
mdl_null = fitglm(tbl, 'y ~ 1', 'Distribution', 'binomial');
LL_full = mdl_full.LogLikelihood;
LL_null = mdl_null.LogLikelihood;
D = -2 * (LL_null - LL_full); % or: 2*(LL_full - LL_null)
df = mdl_full.NumEstimatedCoefficients - mdl_null.NumEstimatedCoefficients;
p_slope_cycle_num = 1 - chi2cdf(D, df);
slope_chi_square = D;
ylim([0 0.5])


%% Are TIPs tighter around alpha-star than non-TIPs?
% Flags, for optional changes
exclude_final_cycle = 0;
consider_data_after_approach_flag = 1;

%Extract alpha values for each group of TIP/cont observations:
trial_num_vec = [res.trial_num];
side_num_vec = [res.side_num];
t_vec = res.tt.contact;

alpha_cycle_vals = cell(0);
is_tip = [];
is_last = [];
pump_num_vec = [];
for trial_idx = 1:length(trial_num_vec)
    trial_num = trial_num_vec(trial_idx);
    side_num = side_num_vec(trial_idx);
    trim_idx = allJobs{trial_num}.Tracks.trim_idx;
    if isnan(trim_idx)
        continue
    end

    d_ang = allJobs{trial_num}.Tracks.diffAngles;
    d_curv = allJobs{trial_num}.Tracks.diffCurv;
    alpha_vals = cart2pol(d_ang,d_curv);

    dd = allJobs{trial_num}.Tracks.whiskerObjDist<cont_thresh_px;
    cont = dd./dd;
    side = mode(allJobs{trial_num}.Tracks.side,2);
    side = side == side_num;
    side = side./side;
    cont = cont.*side;

    cyc_t = t_vec{trial_idx};
    if ~consider_data_after_approach_flag
        cyc_t(end) = allJobs{trial_num}.Tracks.trim_idx;
    end

    for cyc_num = 1:((size(cyc_t,1))-exclude_final_cycle)
        tt = cyc_t(cyc_num,:);
        tt = tt(1):tt(end);

        cyc_alpha = cont(:,tt).*alpha_vals(:,tt);
        cyc_alpha = cyc_alpha(~isnan(cyc_alpha));

        alpha_cycle_vals{end+1} = cyc_alpha;

        pump_idx = find([pump_struct.trial_num] == trial_num & ...
            [pump_struct.side_num] == side_num & ...
            [arrayfun(@(x) x.start_end(1),pump_struct)] == tt(1));
        is_tip(end+1) = pump_struct([pump_idx]).is_tip;
        is_last(end+1) = cyc_num == (size(cyc_t,1));
        pump_num_vec(end+1) = pump_idx;
    end
end
is_tip = logical(is_tip);

% Calculate MAD
% dispersion of alpha around alpha-star
alpha_func = @(x) median(min(...
    abs(circ_dist(x,pi/2-alpha_star_displacement)),...
    abs(circ_dist(x,-pi/2-alpha_star_displacement))));
x_rng = -0.1:0.1:pi/2;
alpha_cycle_pop = cellfun(@(x) alpha_func(x),alpha_cycle_vals);

% Depending on exclusion criteria, contacts may rarely not be identified
% correctly
nan_idx = isnan(alpha_cycle_pop);
alpha_cycle_pop = alpha_cycle_pop(~nan_idx);
is_tip = is_tip(~nan_idx);

% Plot
fh(9) = axes('Units','Centimeters','Position',[6 17 2.5 2.5]);
[y,x] = histcounts(alpha_cycle_pop(~is_tip),x_rng,'Normalization','cdf');
x = x(1:end-1)+0.5*mode(diff(x));
plot(x,y,'r-')
hold on
[y,~] = histcounts(alpha_cycle_pop(is_tip),x_rng,'Normalization','cdf');
plot(x,y,'r-o');
ylim([0 1])
xlim([0 pi/2])
xlabel('$MAD(\alpha;\alpha^{*}))$','Interpreter','latex');
ylabel('CDF')
set(gca,'XTick',linspace(0,pi/2,5),'XTickLabel',makePolarTicks(linspace(0,pi/2,5)),'TickLabelInterpreter','latex')

% Get statistics
[p_tightness_of_tips,~,p_tight_stats] = ranksum(alpha_cycle_pop(~is_tip),alpha_cycle_pop(is_tip),'tail','right');
p_tightness_of_tips = [p_tightness_of_tips,p_tight_stats.ranksum,sum(~is_tip),sum(is_tip)];

% Add legend
leg = legend({'TIP';'Contact'},'Units','centimeters');
leg.Position = [7.5 17.5 2.25 0.8];

%% Within cycle convergence? (using dispersion minimization)
% Make sure you don't exclude anything, start again with data:
trial_num = [res.trial_num];
side_num = [res.side_num];
tt = [res.tt.contact];

trial_num_vec = [pump_struct.trial_num];
side_num_vec = [pump_struct.side_num];
t_vec = [pump_struct.start_end];
t_vec = reshape(t_vec,2,length(t_vec)/2)';

% Where do TIPs occur in the sequence
is_tip_vec = cell(0);
for k = 1:length(trial_num)
    is_tip_temp = [];
    for kk = 1:size(tt{k},1)
        trial_cond = trial_num_vec == trial_num(k);
        side_cond = side_num_vec == side_num(k);
        t_cond = all(t_vec == tt{k}(kk,:),2);
        idx = find(trial_cond & side_cond & t_cond');
        is_tip_temp(end+1) = pump_struct(idx).is_tip;
    end
    is_tip_vec{k} = is_tip_temp;
end

% Collect TIP data
tip_cyc = [];
for kk = 1:length(is_tip_vec)
    for k = 1:(length(is_tip_vec{kk}))
        if is_tip_vec{kk}(k)
            cyc_tt = tt{kk}(k,:);
            tip_cyc(end+1,:) = [res.trial_num(kk),res.side_num(kk),cyc_tt, k == length(is_tip_vec{kk})];
        end
    end
end
if exclude_final_cycle % This flag was set in the previous analysis
    tip_cyc(logical(tip_cyc(:,end)),:) = [];
end


% Extract values for both protractions;
% Same as above, but only for protraction:
tip_pairs = cell(0);
counter = 0;
for k = 1:size(tip_cyc,1)
    trial_num = tip_cyc(k,1);
    side_num = tip_cyc(k,2);
    t = tip_cyc(k,3):tip_cyc(k,4);

    trial_cond = [pump_struct.trial_num] == trial_num;
    side_cond = [pump_struct.side_num] == side_num;
    t_cond = [pump_struct.start_end];
    t_cond = t_cond(1:2:end) == t(1);
    idx = find(trial_cond & side_cond & t_cond);

    pump_idx = pump_struct(idx).trace_signs == 1;
    first_prot_idx = find(pump_idx);
    sec_prot_idx = first_prot_idx(2);
    first_prot_idx = first_prot_idx(1);
    signed_traces = pump_struct(idx).signed_traces;
    first_prot_start = find(~isnan(signed_traces(first_prot_idx,:)),1);
    sec_prot_start = find(~isnan(signed_traces(sec_prot_idx,:)),1);


    pump_idx = pump_struct(idx).trace_signs ~= 1;
    first_ret_idx = find(pump_idx);
    sec_ret_idx = first_ret_idx(2);
    first_ret_idx = first_ret_idx(1);

    first_ret_start = find(~isnan(signed_traces(first_ret_idx,:)),1);
    sec_ret_start = find(~isnan(signed_traces(sec_ret_idx,:)),1);

    % For start of protraction to start of retraction
    t1 = first_prot_start:(first_ret_start-1);
    t2 = sec_prot_start:(sec_ret_start-1);

    t_vals = {t1; t2};
    tip_pair = cell(0);
    bad_flag = false;
    for mm = 1:length(t_vals)
        t_val = t_vals{mm};
        d_ang = allJobs{trial_num}.Tracks.diffAngles(:,t_val);
        d_curv = allJobs{trial_num}.Tracks.diffCurv(:,t_val);
        cont = allJobs{trial_num}.Tracks.whiskerObjDist(:,t_val)<11.507;
        cont = cont./cont;
        side = mode(allJobs{trial_num}.Tracks.side,2) == side_num;
        side = side./side;
        alpha_vals = cart2pol(d_ang,d_curv);
        alpha_vals = alpha_vals.*cont.*side;
        alpha_vals = alpha_vals(~isnan(alpha_vals));
        if isempty(alpha_vals)
            counter = counter+1;
            bad_flag = true;
        end
        tip_pair{mm} = alpha_vals;
    end
    if ~bad_flag
        tip_pairs{end+1} = tip_pair;
    end
end

% Get representative alpha-value
first_pump = cellfun(@(x) dispersionMinimization(x{1}),tip_pairs);
sec_pump = cellfun(@(x) dispersionMinimization(x{2}),tip_pairs);

% Plot:
fh(16) = axes('Units','Centimeters','Position',[11.5 17 2.5 2.5]);
alpha_star_val = pi/2-alpha_star_displacement;
[y,x] = histcounts(abs(first_pump-alpha_star_val),linspace(0,pi,101),'Normalization','cdf');
x = x(1:end-1)+0.5*mode(diff(x));
plot(x,y,'Color',[0 0 0])
hold on
[y,x] = histcounts(abs(sec_pump-alpha_star_val),linspace(0,pi,101),'Normalization','cdf');
x = x(1:end-1)+0.5*mode(diff(x));
plot(x,y,'Color',[0 0 1])
x_tick = sort([setdiff(linspace(0,pi,5),pi/2),pi/2-alpha_star_displacement]);
x_tick_lbl = makePolarTicks(linspace(0,pi,5));
set(gca,'XTick',x_tick,'XTickLabel',x_tick_lbl,'TickLabelInterpreter','latex')
hold on
ylabel('CDF');
xlabel('$|\arg\min \mathrm{MAD}(\alpha)-\alpha^{*}|$','Interpreter','latex')
leg = legend({'1st protraction';'2nd protraction'},'Units','centimeters');
leg.Position = [13.25 17.5 3 1];
xlim([0 1.7]);

% Statistics:
[within_cyc_stats(1),~,pumps_wilcoxson_stats] = signrank(abs(first_pump-(pi/2-alpha_star_displacement)),abs(sec_pump-(pi/2-alpha_star_displacement)),'tail','right');
within_cyc_stats(2) = pumps_wilcoxson_stats.signedrank;
within_cyc_stats(3:4) = [length(first_pump),length(sec_pump)];


%% TIP-gating analysis:

% Only take data during protraction period of the contact
just_prot_flag = true; 

% Prepare to extract cycle data
trial_num = [res.trial_num];
side_num = [res.side_num];
tt = [res.tt.contact];
preceding_alpha_vals = cell(0);
self_is_tip = [];
next_is_tip = [];
just_prot_idx = cell(0);

for trial_idx = 1:length(trial_num)
    d_ang = allJobs{trial_num(trial_idx)}.Tracks.diffAngles;
    d_curv = allJobs{trial_num(trial_idx)}.Tracks.diffCurv;
    cont = allJobs{trial_num(trial_idx)}.Tracks.whiskerObjDist<cont_thresh_px;
    cont = cont./cont;
    side = mode(allJobs{trial_num(trial_idx)}.Tracks.side,2) == side_num(trial_idx);
    side = side./side;
    d_ang = d_ang.*side.*cont;


    pump_trial = [pump_struct.trial_num] == trial_num(trial_idx);
    pump_side = [pump_struct.side_num] == side_num(trial_idx);
    pump_t = reshape([pump_struct.start_end],2,length([pump_struct.start_end])/2)';


    [alpha_vals,r_vals] = cart2pol(d_ang,d_curv);
    for kk = 1:(length(is_tip_vec{trial_idx})-1)

        cyc_t = allJobs{trial_num(trial_idx)}.Tracks.cycles{side_num(trial_idx)}.cycle_start_end;
        cyc_idx = find(tt{trial_idx}(kk,1) == cyc_t(:,1));
        if isempty(cyc_idx)
            disp('Bad frame match?')
            continue
        end
        cyc_alpha = alpha_vals(:,cyc_t(cyc_idx,1):cyc_t(cyc_idx,2));
        cyc_r = r_vals(:,cyc_t(cyc_idx,1):cyc_t(cyc_idx,2));
        pump_idx = find(pump_t(:,1)' == cyc_t(cyc_idx,1) & pump_trial & pump_side);
        % Mark when protraction occurs:
        prot_t = (pump_struct(pump_idx).trace_signs == 1).*(pump_struct(pump_idx).signed_traces);
        prot_t = prot_t./prot_t;
        [~,col] = find(~isnan(prot_t));

        % Mark which values are from protraction:
        prot_t = prot_t(:,cyc_t(cyc_idx,1):cyc_t(cyc_idx,2));
        prot_t = any(~isnan(prot_t),1);
        prot_t = prot_t.*ones(size(cyc_alpha));
        prot_t = prot_t(~isnan(cyc_alpha));

        cyc_alpha = cyc_alpha(~isnan(cyc_alpha));
        cyc_r = cyc_r(~isnan(cyc_r));

        if isempty(cyc_alpha)
            continue
        end
        preceding_alpha_vals{end+1} = cyc_alpha;
        just_prot_idx{end+1} = prot_t; 

        self_is_tip(end+1) = is_tip_vec{trial_idx}(kk);
        next_is_tip(end+1) = 0;
        if is_tip_vec{trial_idx}(kk+1)
            next_is_tip(end) = 1;
        end

    end
end

% Get the two types of  cycle labels:
before_cont = ~self_is_tip & ~next_is_tip; % Cont followed by cont
before_tip = ~self_is_tip & next_is_tip; % Cont followed by TIP
conds = {[before_cont;before_tip]};

% Only protraction periods:
if just_prot_flag
    preceding_alpha_vals = cellfun(@(x,y) x(logical(y)),preceding_alpha_vals,just_prot_idx,'Uni',0);
end

% Function to assign representative values:
func = @(x) abs(dispersionMinimization(x)-(pi/2-alpha_star_displacement));

% Remove empty cycles:
rmv_idx = cellfun(@(x) length(x) < 1,preceding_alpha_vals);
preceding_alpha_vals(rmv_idx) = [];
conds = cellfun(@(x) x(:,~rmv_idx),conds,'Uni',0);
conds = conds{1};

% Alpha values, first before cont, followed by before TIP
cyc_vals = [preceding_alpha_vals(conds(1,:) == 1),...
    preceding_alpha_vals(conds(2,:) == 1)];
% Apropriate labels
cyc_labels = [zeros(1,sum(conds(1,:) == 1)),...
    ones(1,sum(conds(2,:) == 1))];
% Get representative values
cyc_vals = cellfun(@(x) func(x) ,cyc_vals);

% Sort values and labels, according to the values
[~,sort_idx] = sort(cyc_vals);
cyc_vals = cyc_vals(sort_idx);
cyc_labels = cyc_labels(sort_idx);
tip_gating_Ns = [sum(cyc_labels),sum(~cyc_labels)]; % Ns, contacts followed by: TIPs, contacts

% Plot probability, using a moving average window:
axes('Units','centimeters','Position',[1.5 17 2.5 2.5])
NN = 30;
xx = movmean(cyc_vals,NN,'Endpoints','discard');
yy = movmean(cyc_labels,NN,'Endpoints','discard');
pp1 = plot(xx,yy,'.','Color',[1 1 1]*0.5);
ylim([0 0.7])
hold on
set(gca,'XTick',linspace(0,pi/2,5),'XTickLabel',makePolarTicks(linspace(0,pi/2,5)),'TickLabelInterpreter','latex')
ylabel('p(next cycle is TIP)')
xlabel('$|\arg\min \mathrm{MAD}(\alpha)-\alpha^{*}|$','Interpreter','latex');
drawnow

% Fit exponential model:
model = @(p, x) p(1) * exp(p(2) * x) + p(3);
x0 = [1, -1, 0];  % Initial guess
[beta_fit, ~, ~, ~, ~, ~, ~] = ...
    lsqcurvefit(model, x0, xx, yy);

% Get adjusted r-squared:
n = numel(yy);     % number of data points
k = numel(beta_fit);       % number of parameters
y_pred = model(beta_fit, xx);      % predicted values
residuals = yy - y_pred;
SS_res = sum(residuals.^2);                      % Residual sum of squares
SS_tot = sum((yy - mean(yy)).^2);        % Total sum of squares
R2 = 1 - SS_res / SS_tot;                        % R-squared
adjR2 = 1 - (1 - R2)*(n - 1)/(n - k);            % Adjusted R-squared

% Get p-value vs constant
df1 = k - 1;
df2 = n - k;
MS_model = (SS_tot - SS_res) / df1;
MS_resid = SS_res / df2;
F = MS_model / MS_resid;
exp_p_value = 1 - fcdf(F, df1, df2);
pp2 = plot(xx,model(beta_fit,xx),'k--','LineWidth',1.5);


%% Tidy up:
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
%%
% A bit more tidying
axesHandles = findall(gcf, 'Type', 'axes');
legendHandles = findall(gcf, 'Type', 'legend');

translate_down_x_cm = 8;

% Move axes
for ax = axesHandles'
    ax.Units = 'centimeters';
    ax.Position(2) = ax.Position(2) - translate_down_x_cm;
end

% Move legends
for lgd = legendHandles'
    lgd.Units = 'centimeters';
    lgd.Position(2) = lgd.Position(2) - translate_down_x_cm;
end
drawnow expose

%% Fig 5 stats
disp('Fig 5 stats:')

disp(' ');
disp('Fig 5b+2c:')
helper_ttls = {'Single-pump free-air';'Contact';'Double-pump';'TIP'};
for k = 1:4
    disp([helper_ttls{k},': n = ',num2str(length(all_cyc_duration{k}))])
end
disp(' ');
disp('TIP duration vs cont, one-tail, Wilcoxon''s rank-sum:')
disp(['n_cont = ',num2str(length(all_cyc_duration{3})),'n_tip = ',num2str(length(all_cyc_duration{4})),';']);
disp([' p = ',num2str(p_val_tip_vs_cont_dur),'; W = ',num2str(stats_tip_vs_cont_dur.ranksum)]);
disp(['TIP is longer-duration by: ',num2str(TIP_elongated_prct-1)]);

disp(' ');
disp(['N_trials for rest = ',num2str(length([res.trial_num]))]);

disp(' ');
disp('Fig 5d:')
disp('Logistic regression vs null (constant):');
disp(['Chi-square: ',num2str(slope_chi_square)]);
disp([' p = ',num2str(p_slope_cycle_num)]);
disp(['Last cycle vs previous, one-tail Fisher''s exact test:']);
disp(['non_last_N = ',num2str(length(grp2)),'; last_N = ',num2str(length(grp1)),';']);
disp(['non_last_p = ',num2str(mean(grp2)),'; last_p = ',num2str(mean(grp1)),';']);
disp(['p = ',num2str(p_binomial_last_cycle(1))]);
disp('Ns from last backwards:')
for k = 1:length(tip_bars_Ns)
    disp(num2str(tip_bars_Ns(k)))
end


disp(' ');
disp('Fig 5e, TIP emergence, exp. fit:')
disp(['N_cycles = ',num2str(length(cyc_vals)),'; mov-avg_N = 30;']);
disp('N cont followed by TIP:');
disp(num2str(tip_gating_Ns(1)));
disp('N cont followed by cont:');
disp(num2str(tip_gating_Ns(2)));
disp(['adj_R^2 = ',num2str(adjR2),'; p = ',num2str(exp_p_value)]);


disp(' ');
disp('Fig 5f:')
disp('Tightness of TIP vs cont, one-tail Wilcoxon''s rank-sum:')
disp(['n_TIP = ',num2str(p_tightness_of_tips(4)),'; n_Cont = ',num2str(p_tightness_of_tips(3)),';']);
disp([' p = ',num2str(p_tightness_of_tips(1)),'; W = ',num2str(p_tightness_of_tips(2))]);


disp(' ');
disp('Fig 5g:')
disp('First vs second pump:')
disp(['n_first = ',num2str(within_cyc_stats(3)),'; n_sec = ',num2str(within_cyc_stats(4)),';']);
disp('one-tail Wilcoxon signed-rank(paired), deviations from alpha-star:')
disp([' p = ',num2str(within_cyc_stats(1)),'; W = ',num2str(within_cyc_stats(2))]);


end