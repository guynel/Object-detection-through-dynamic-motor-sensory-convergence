function [] = generateSupp3(allJobs,px2mm,extend_to_cycle_end_flag,cont_thresh_px,alpha_star_val)
%% Plot extended Fig. 3

%% Get cycle data:
side_pref_thresh = 0.9;
pre_cont_thresh_coeff = 1;

res = getDataPerCycle(allJobs,px2mm,side_pref_thresh,pre_cont_thresh_coeff,cont_thresh_px,extend_to_cycle_end_flag);

%% 
output = [];
for k = 1:(length(res.contact)) % Iterate over trials
    if size(res.tt.contact{k},1) < 2 % Only it has 3 or more cycles, so we can do lagged correlation without including last cycle
        continue
    end
    side_num = res.side_num(k);
    trial_num = res.trial_num(k);
    % Head-center, taken with respect to object center
    hc = allJobs{trial_num}.Tracks.headCenter;
    sc = allJobs{trial_num}.stimContour.stimCenter;
    hc = hc-sc';

    for kk = 1:(size(res.tt.contact{k},1)-1) % For each cycle and the following, excluding the pairs of cycles that are penultimate-to-last 
        tt_curr = res.tt.contact{k}(kk,:);
        tt = allJobs{trial_num}.Tracks.cycles{side_num}.cycle_start_end;
        cyc_idx = find(all(tt == tt_curr,2));
        tt_next = tt(cyc_idx+1,:);

        % Change in the median head-position between cycles
        tt_curr = tt_curr(1):tt_curr(end);
        tt_next = tt_next(1):tt_next(end);
        hc_curr = median(hc(:,tt_curr)');
        hc_next = median(hc(:,tt_next)');
        d_hc = (hc_next)-(hc_curr);

        % Pair with alpha value in the first of the cycles:
        alpha_val = dispersionMinimization(res.contact{k}{kk});
        output = [output; [alpha_val,d_hc(:,2)]];
    end
end

% Plot, use moving window to smooth sorted observations:
fig = makeFullPagePDF();
q = [abs(circ_dist(output(:,1),alpha_star_val)),0.1*(1/px2mm)*abs(output(:,2))];
[~,idx] = sort(q(:,1));
q = q(idx,:);
sc = scatter(q(:,1),q(:,2),10,[0 0 1],'filled', 'MarkerEdgeAlpha',0.5,'MarkerFaceAlpha',0.5);
ylabel({'$|\Delta median(head-center_{y})|$';'(from cycle n to n+1) $[cm]$'},'Interpreter','latex');
xlabel('dist$(\alpha(n),\alpha^{*})$','Interpreter','latex');
xlim([0 pi-alpha_star_val])
hold on
NN = 30;
p_red = plot(movmean(q(:,1),NN),movmean(q(:,2),NN),'r','LineWidth',1.5);
sim_res = [];
for iter = 1:1000
    qq = q;
    sim_res(:,end+1) = movmean(q(randperm(length(q),length(q)),2),NN);
end
pp = patch([q(:,1);flipud(q(:,1))],[prctile(sim_res,2.5,2);...
    prctile(sim_res,97.5,2)],[0 0 0],'FaceAlpha',0.1);
uistack(pp,'bottom');
set(gca,'FontSize',12)
set(gca,'Units','centimeters','Position',[4 10 4 4])
leg = legend([sc,p_red],{['Cycles (n = ',num2str(length(q)) ')']; 'Moving-average (window_n = 30)'});
leg.Location = 'northoutside';
leg.ItemTokenSize = [10,1];
leg.Units = 'centimeters';
leg.Position = [5 14.5 3 1];
set(gca,'XTick', linspace(0,pi/2,5),'XTickLabel',makePolarTicks(linspace(0,pi/2,5)),'TickLabelInterpreter','latex')

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