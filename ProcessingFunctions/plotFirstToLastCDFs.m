function [p_vals,p_handles] = plotFirstToLastCDFs(plot_var,plot_func,panel_holder,x_rng,x_lbls,first_and_last_flag)
% Used to plot the top row of Fig. 4

% Plotting CDFs
clrs = [0.90, 0.62, 0; 0 0 0; 0.94, 0.89, 0.26; 0.72, 0.53, 0.04];
clrs = [1, 0.65, 0; 0 0 0; 1, 0, 0; 0, 1, 0];

var_data = plot_var;
% Single-contact cycle data
first_and_last = var_data(cellfun(@(x) length(x)==1,var_data));
first_and_last = [first_and_last{:}];

var_data = plot_var;
var_data = var_data(cellfun(@(x) length(x)>1,var_data));
first = cellfun(@(x) x{1},var_data,'Uni',0); % First
last = cellfun(@(x) x{end},var_data,'Uni',0); % Last
interim = var_data(cellfun(@(x) length(x)>2,var_data));
interim = cellfun(@(x) x(2:(end-1)),interim,'Uni',0);
interim = [interim{:}]; % Interim

% Get representative values
first_and_last_vals = cellfun(@(x) plot_func(x),first_and_last);
first_vals = cellfun(@(x) plot_func(x),first);
interim_vals = cellfun(@(x) plot_func(x),interim);
last_vals = cellfun(@(x) plot_func(x) ,last);

var_rng = x_rng; % Plotting range

% Plot:
hold on
[y_first,x] = histcounts(first_vals,var_rng,'Normalization','cdf');
x = x(1:end-1)+0.5*mode(diff(x));
plot(x,y_first,'Color',clrs(1,:));
if first_and_last_flag
    [y_interim,~] = histcounts(interim_vals,var_rng,'Normalization','cdf');
    plot(x,y_interim,'Color',clrs(2,:));
end
[y_last,~] = histcounts(last_vals,var_rng,'Normalization','cdf');
plot(x,y_last,'Color',clrs(3,:));
[y_single,~] = histcounts(first_and_last_vals,var_rng,'Normalization','cdf');
plot(x,y_single,'Color',clrs(4,:));
xlim([var_rng(1),var_rng(end)]);
ylim([0 1]);
ylabel('CDF');
xlabel(x_lbls,'Interpreter','latex');

rng(121987,"twister"); % Set random seed for reproducibility of bootstrap

% Permutation test for differenece between single and last
sim_res = [];
for iter = 1:1000
    grpd = [last_vals,first_and_last_vals];
    shuff_idx = randperm(length(grpd));
    grp1 = grpd(shuff_idx(1:length(last_vals)));
    grp2 = grpd(shuff_idx((length(last_vals)+1):end));
    sim_res(iter) = abs(median(grp1)-median(grp2));
end
gt = abs(median(last_vals)-median(first_and_last_vals));
p2 = min(mean(sim_res<gt),mean(sim_res>gt));

% Permutation test for differenece between single and first
sim_res = [];
for iter = 1:1000
    grpd = [first_vals,first_and_last_vals];
    shuff_idx = randperm(length(grpd));
    grp1 = grpd(shuff_idx(1:length(first_vals)));
    grp2 = grpd(shuff_idx((length(first_vals)+1):end));
    sim_res(iter) = abs(median(grp1)-median(grp2));
end
gt = abs(median(first_vals)-median(first_and_last_vals));
p1 = min(mean(sim_res<gt),mean(sim_res>gt));

% Plotting parameters for significance bars
p_vals = [p1,p2];
x_rng = xlim();
y_rng = ylim();
p_handles_x_start = (range(x_rng)*[0.7])+min(x_rng);
p_handles_x_width = (range(x_rng)*[0.075]);
p_handles_y_start = (range(y_rng)*[0.075])+min(y_rng);
p_handles_y_width = (range(y_rng)*[0.05]);

% Plotting parameters for significance bars: single vs first
p1_handles(1) = rectangle('Position',[p_handles_x_start,p_handles_y_start,...
    p_handles_x_width,p_handles_y_width],'FaceColor',clrs(4,:),'EdgeColor','k');
p1_handles(2) = rectangle('Position',[p_handles_x_start+p_handles_x_width,p_handles_y_start,...
    p_handles_x_width,p_handles_y_width],'FaceColor',clrs(1,:),'EdgeColor','k');
p1_handles(3) = text(p_handles_x_start+2.75*p_handles_x_width,...
    p_handles_y_start*0.9,'*','FontSize',16,'Color','k',...
    'HorizontalAlignment','center','VerticalAlignment','middle');


% Plotting parameters for significance bars: single vs last
p2_handles(1) = rectangle('Position',[p_handles_x_start,p_handles_y_start+4*p_handles_y_width,...
    p_handles_x_width,p_handles_y_width],'FaceColor',clrs(4,:),'EdgeColor','k');
p2_handles(2) = rectangle('Position',[p_handles_x_start+p_handles_x_width,p_handles_y_start+4*p_handles_y_width,...
    p_handles_x_width,p_handles_y_width],'FaceColor',clrs(3,:),'EdgeColor','k');
p2_handles(3) = text(p_handles_x_start+2.75*p_handles_x_width,...
    p_handles_y_start*0.9+4*p_handles_y_width,'*','FontSize',16,'Color','k',...
    'HorizontalAlignment','center','VerticalAlignment','middle');
p_handles = {p1_handles;p2_handles};

% We output handles and p-value, perform FDR correction, and delete bars
% that do not survive correction

end
