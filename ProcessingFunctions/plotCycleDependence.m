function [p_vals,lin_reg_plot] = plotCycleDependence(plot_var,var_rng,plot_func,mult_fact,varargin)
% Used to plot the bottom row of Fig. 4
hold on

% Get first and last cycle data for each relevant trial
all_first = cellfun(@(x) x{1},plot_var,'Uni',0);
all_first_vals = cellfun(@(x) plot_func(x),all_first);
all_last = cellfun(@(x) x{end},plot_var,'Uni',0);
all_last_vals = cellfun(@(x) plot_func(x),all_last);
N_cycles = cellfun(@(x) length(x),plot_var);
N_cycles(N_cycles>5) = 5; % If more than 5 contact occur, lump it into the >= 5 group

% Fit a linear regression for alpha value in first cycle, according to
% total number of cycles in trial, and plot it
mdl = fitlm(N_cycles', all_first_vals');
p_values = mdl.Coefficients.pValue;
p_vals(1) = p_values(2);
clrs = [1.00, 0.65, 0.00; 1, 0, 0];
lin_reg_plot(1) = plot(0:0.01:6,predict(mdl,(0:0.01:6)'),'--','LineWidth',2,'Color',[1 1 1]*0.5);

% We output the p-value and regression-plot handles; then after FDR
% correction we remove those the insignificant

% Fit a linear regression for alpha value in last cycle, according to
% total number of cycles in trial, and plot it
mdl = fitlm(N_cycles', all_last_vals');
p_values = mdl.Coefficients.pValue;
p_vals(2) = p_values(2);
lin_reg_plot(2) = plot(0:1:6,predict(mdl,(0:1:6)'),'--','LineWidth',2,'Color',[1 1 1]*0);
xlim([0 6])
set(gca,'XTick',1:6);

% Plot underlying data as well
sc = gobjects(0);
for cycle_num = 1:max(N_cycles)
    % First cycles:
    cycle_obs = all_first_vals(N_cycles == cycle_num);
    [y,x] = ksdensity(cycle_obs,var_rng,'Function','pdf',varargin{1}{:});
    x = [x,x(end),x(1),x(1)];
    y = [y,0,0,y(1)];
    y = y./sum(y);
    y = y*mult_fact;
    pp = patch(cycle_num+y,x,clrs(1,:),'FaceAlpha',[0.4],'LineWidth',1);
    uistack(pp,'bottom')

    % For cycle num = 1, first and last will be same data, so when the
    % mean is scatterplotted on the data, one marker will cover the other;
    % here we account for that by having one marker within the other
    if cycle_num == 1
        sc_top = scatter(cycle_num,mean(cycle_obs),8,'o','filled','MarkerFaceColor',clrs(1,:));
    else
        sc(end+1) = scatter(cycle_num,mean(cycle_obs),'o','filled','MarkerFaceColor',clrs(1,:),'MarkerEdgeColor','k','LineWidth',1);
    end
    
    % Last cycles
    cycle_obs = all_last_vals(N_cycles == cycle_num);
    [y,x] = ksdensity(cycle_obs,var_rng,'Function','pdf',varargin{1}{:});
    x = [x,x(end),x(1),x(1)];
    y = [y,0,0,y(1)];
    y = y./sum(y);
    y = y*mult_fact;
    pp = patch(cycle_num-y,x,clrs(2,:),'FaceAlpha',[0.4],'LineWidth',1);
    sc(end+1) = scatter(cycle_num,mean(cycle_obs),'o','filled','MarkerFaceColor',clrs(2,:),'MarkerEdgeColor','k','LineWidth',1);
   
    
end
set(gca,'XDir','reverse');

% Order elements to look nicer
y_lim = var_rng;
y_lim = [min(y_lim),max(y_lim)];
ylim(y_lim);
for k = 1:length(lin_reg_plot)
    uistack(lin_reg_plot(k),'top');
end
for k = 1:length(sc)
    uistack(sc(k),'top')
end
uistack(sc_top,'top');
end
