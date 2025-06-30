function [] = generateFig1TraceExample(allJobs,px2mm,fig1_handle)

%% Variables in isolation and in pairs

fontSize = 12;
lineW = 1;
cont_thresh = 11.8509;

% Example:
plot_trial = 54; % Trial
cont_tracks = 2; % Contacting whisker
non_cont_tracks = 1; % Non-contacting whisker
t = 225:571; % Relevant frames

jobTracks = allJobs{plot_trial}.Tracks; % Data for selected job

% Streamline data and labels into a struct so plotting is more
% straightforward
var_names = {'smoothAngles','smoothCurv';'smoothDiffAngles','smoothDiffCurv'}'; % Names of variables within the data struct
var_legend = {'$\theta$ $[degrees]$','$\kappa$ $[m^{-1}]$';'$\dot{\theta}$ $[degrees/s]$','$\dot{\kappa}$ $[m^{-1}/s]$'}'; % Names for legend
res = struct();
for i = 1:length(var_names(:))
    
    
    res.cont_track.(var_names{i}) = nanmedian(jobTracks.(var_names{i})(cont_tracks,:),1); % Contact track
    mask = any(jobTracks.whiskerObjDist(cont_tracks,:)<cont_thresh,1); % Cont frames
    res.cont_track.mask = mask./mask;
    res.cont_track.clrs = [150 150 150; 255 0 0]./255; % Colors for plotting the contact track
    res.non_cont_track.(var_names{i}) = jobTracks.(var_names{i})(non_cont_tracks,:); % Non-contact track
    mask = any(jobTracks.whiskerObjDist(non_cont_tracks,:)<cont_thresh,1); % Cont frame (degenerate, will not pick up anything)
    res.non_cont_track.mask = mask./mask;
    res.non_cont_track.clrs = [0,0,0; 191 191 191]./255; % Colors (only one of these will be used if there are no contact frames for the no contact track)
end
% Data could have NaN values in the edges, which is not nice visually, so
% let's remove those values automatically.
% Start by collecting all relevant data for both tracks in a single matrix
all_plot_data = [];
cont_types = fieldnames(res);
for k = 1:length(cont_types)
    for m = 1:length(var_names(:))
        all_plot_data = [all_plot_data; res.(cont_types{k}).(var_names{m})];
    end
end
% Trim edges:
trim_nan_edges = all(~isnan(all_plot_data));
t_start = find(trim_nan_edges,1);
t_end = find(trim_nan_edges,1,'last');


tt = (t-t(1))*2; % From frames to ms
fig = fig1_handle;
% Plot the individual time series for all four variabls
for m = 1:length(var_names(:))
    
    for k = 1:length(cont_types)
        if k == 1
            fh(m) = axes();
            fh(m).Units = 'Centimeters';
            fh(m).InnerPosition = [2.5 15.25-3.25*(m-1)-0.5*(m>2) 7 1.75];
        end
        var = (1+499*(m>2))*res.(cont_types{k}).(var_names{m})(:,t);
        clr = res.(cont_types{k}).clrs(1,:);
        plot(tt,var,'Color',clr,'LineWidth',lineW);
        mask = res.(cont_types{k}).mask(:,t);
        clr = res.(cont_types{k}).clrs(2,:);
        hold on
        plot(tt,var.*mask,'Color',clr,'LineWidth',lineW);
        if m == length(var_names(:))
            xlabel('Time [ms]','FontSize',12);
        end
        ylabel(var_legend{m},'Interpreter','latex','FontSize',12);
        xlim([tt(1),tt(end)]);
        if m>2 & k >1
            set(gca,'YTickLabels',get(gca,'YTick')./1000);
            curr_ax = gca;
            curr_ax.YLabel.String = {curr_ax.YLabel.String;'x $10^4$'};
        end
        
    end
    if m ~= length(var_names(:))
        continue
    end
    
    y_lim = ylim();
    if m == 1
        ylim([45 145]);
    else
        ylim(1.1*[-1 1].*max(abs(y_lim)));
    end
    
end
% Tidy up
for k = 1:length(fh)
    axes(fh(k));
    line_specs = findall(gca,'Type','Line');
    uistack(line_specs(2),'bottom');
end
axes(fh(1));
line_specs = findall(gca,'Type','Line');
leg = legend(line_specs([2,1,4]),{'Contact'; 'Peri-contact';'Non-contacting whisker'},'Layout','');
leg.Units = 'Centimeters';
leg.FontSize = 9;
leg.FontWeight = 'bold';
leg.Position = [2.5 17 5 1.5];


% Plot the theta-kappa phase-plane
var_names = var_names';
var_legend = var_legend';
max_ax = 0;
k = 1;
for m = length(cont_types):-1:1
    if m == length(cont_types)
        ffh(k) = axes();
        ffh(k).Units = 'Centimeters';
        ffh(k).InnerPosition = [12 12-5*(k-1) 5 5];
    end
    var1 = (1+499*(k>1))*res.(cont_types{m}).(var_names{k,1})(:,t);
    var2 = (1+499*(k>1))*res.(cont_types{m}).(var_names{k,2})(:,t);
    if k == 1
        max_ax_candidate = max(abs([var1(:)-90;var2(:)]));
    else
        max_ax_candidate = max(abs([var1(:);var2(:)]));
    end
    max_ax = max(max_ax,max_ax_candidate);
    clr = res.(cont_types{m}).clrs(1,:);
    plot(var1,var2,'Color',clr,'LineWidth',lineW);
    hold on
    mask = res.(cont_types{m}).mask;
    var1 = var1.*mask(:,t);
    var2 = var2.*mask(:,t);
    clr = res.(cont_types{m}).clrs(2,:);
    plot(var1,var2,'Color',clr,'LineWidth',lineW);

    xlabel(var_legend{k,1},'Interpreter','latex','FontSize',12);
    ylabel(var_legend{k,2},'Interpreter','latex','FontSize',12);

end
% Axis tidy-up
axis equal
if k == 1
    ax = max_ax*[-1 1 -1 1]*1.1+[90 90 0 0];
else
    ax = max_ax*[-1 1 -1 1]*1.1;
end

%%
% Non-contacting whisker during contact - empty panel:
axes('Units','Centimeters','Position',[14.5 5 2.5 2.5])
set(gca,'XTick',[],'YTick',[]);
hold on
plot([-1,1],[0 0],'k--');
plot([0 0],[-1,1],'k--');
axis([-1 1 -1 1])
box on;
% Non-contacting whisker during non-contact - black:
black_line_panel = axes('Units','Centimeters','Position',[12 5 2.5 2.5]);
non_cont_d_theta = 500*res.non_cont_track.smoothDiffAngles(:,t);
non_cont_d_kappa = 500*res.non_cont_track.smoothDiffCurv(:,t);
clr1 = res.non_cont_track.clrs(1,:);
xy1(:,1) = non_cont_d_theta;
xy1(:,2) = non_cont_d_kappa;
plot(xy1(:,1),xy1(:,2),'Color',clr1);
axis([-0.5 0.5 -1 1]*4000);
xlabel([var_legend{2},' x $10^4$'],'Interpreter','latex');
ylabel([var_legend{4},' x $10^4$'],'Interpreter','latex');
hold on
plot(xlim,[0 0],'k--');
plot([0 0],ylim,'k--');
set(black_line_panel,'XTickLabels',get(black_line_panel,'XTick')./1000);
set(black_line_panel,'YTickLabels',get(black_line_panel,'YTick')./1000);
text(-0.65, 0.5, {'Non-contacting';'whisker'}, 'Rotation', 90, ...
    'VerticalAlignment', 'middle', 'HorizontalAlignment', 'center', 'Units', 'normalized', 'FontSize', 10);
% Contacting whisker during non-contact - grey:
axes('Units','Centimeters','Position',[12 7.5 2.5 2.5])
cont_d_theta = 500*res.cont_track.smoothDiffAngles(:,t);
cont_d_kappa = 500*res.cont_track.smoothDiffCurv(:,t);
cont_mask = res.cont_track.mask(:,t);
cont_mask = isnan(cont_mask)./isnan(cont_mask);
clr2 = res.cont_track.clrs(1,:);
xy2(:,1) = cont_d_theta.*cont_mask;
xy2(:,2) = cont_d_kappa.*cont_mask;
plot(xy2(:,1),xy2(:,2),'Color',clr2);
axis([-0.5 0.5 -1 1]*4000);
set(gca,'XTick',[],'YTick',[]);
hold on
plot(xlim,[0 0],'k--');
plot([0 0],ylim,'k--');
title('Non-contact','FontSize',10)
text(-0.65, 0.5, {'Contacting';'whisker'}, 'Rotation', 90, ...
    'VerticalAlignment', 'middle', 'HorizontalAlignment', 'center', 'Units', 'normalized', 'FontSize', 10);
% Contacting whisker during contact - red:
axes('Units','Centimeters','Position',[14.5 7.5 2.5 2.5])
cont_d_theta = 500*res.cont_track.smoothDiffAngles(:,t);
cont_d_kappa = 500*res.cont_track.smoothDiffCurv(:,t);
cont_mask = res.cont_track.mask(:,t);
clr3 = res.cont_track.clrs(2,:);
xy3(:,1) = cont_d_theta.*cont_mask;
xy3(:,2) = cont_d_kappa.*cont_mask;
plot(xy3(:,1),xy3(:,2),'Color',clr3);
set(gca,'XTick',[],'YTick',[]);
axis([-0.5 0.5 -1 1]*4000);
hold on
plot(xlim,[0 0],'k--');
plot([0 0],ylim,'k--');
title('Contact','FontSize',10);
% Inset: MAD-minimized ellipses for each sub-panel, overlaid:
axes('Units','Centimeters','Position',[16.6 7 2 2])
var_pairs = {xy1;xy2;xy3};
var_clrs = [clr1; clr2; clr3];
for k = 1:length(var_pairs)
    nan_idx = any(isnan(var_pairs{k}),2);
    plot_mad_aligned_ellipse(var_pairs{k}(~nan_idx,1),var_pairs{k}(~nan_idx,2),1,100,var_clrs(k,:))
    hold on
end
% Tidy up:
set(gca,'XTick',[],'YTick',[]);
axis([-0.5 0.5 -1 1]*2);
plot(xlim,[0 0],'k--');
plot([0 0],ylim,'k--');
scatter(0,0,'k.')

%% Inset: illustrate polar coordinates:

axes('Units','Centimeters','Position',[16.6 3.5 2 2])
hold on
plot([-1 1],[0,0],'k','LineWidth',1)
plot([0,0],[-1 1],'k','LineWidth',1)
axis([-1 1 -1 1]*1.1)
axis equal
set(gca,'XTick',[],'YTick',[])

angle = 40;
r = 1;
theta = linspace(0, deg2rad(angle), 100);
x = [0, r * cos(theta)];
y = [0, r * sin(theta)];
fill(x, y, 'k', 'FaceAlpha', 0.1)  % Blue color with transparency
quiver(0,0,1.2*r*cos(deg2rad(angle)), 1.2*r*sin(deg2rad(angle)),'Color',[0.5 0.5 0.5],'LineWidth',2,'MaxHeadSize',10)
text(0.5, 0.3, '\alpha', 'FontSize', 10,'FontWeight','bold','Color',[0.5 0.5 0.5]);
text(0.9,0.9, 'r', 'FontSize', 10,'FontWeight','bold','Color',[0.5 0.5 0.5]);
axis([-1,1,-1,1]*1.2)

%% General tidy-up:

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



end