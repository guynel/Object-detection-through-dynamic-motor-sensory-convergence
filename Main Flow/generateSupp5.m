function [] = generateSupp5(allJobs)
%% Extended Fig. 5

%% Plot example
trial_num = 95; % Example Trial

% Set up page:
makeFullPagePDF();
axes('Units','Centimeters','Position',[1.5 20 12 3])
ang = allJobs{trial_num}.Tracks.origAngles;
pos_ang = ang.*((ang>0)./(ang>0));
tt = repmat(1:size(pos_ang,2),size(pos_ang,1),1);
scatter(tt(~isnan(pos_ang)),pos_ang(~isnan(pos_ang)),3,'k','filled','MarkerFaceAlpha',1,'MarkerEdgeAlpha',1)
set(gca,'XTickLabel','');
ylabel('\theta_{right}')

% Get angles and scatter them as first panel
axes('Units','Centimeters','Position',[1.5 16.5 12 3])
ang = allJobs{trial_num}.Tracks.origAngles;
neg_ang = ang.*((ang<0)./(ang<0));
tt = repmat(1:size(neg_ang,2),size(neg_ang,1),1);
scatter(tt(~isnan(neg_ang)),abs(neg_ang(~isnan(neg_ang))),3,'k','filled','MarkerFaceAlpha',1,'MarkerEdgeAlpha',1)
xlabel('Time [ms]')
ylabel('\theta_{left}')

% Repeat, but plot smoothed angle traces on them as well
% For positive angles:
axes('Units','Centimeters','Position',[1.5 11 12 3])
ang = allJobs{trial_num}.Tracks.origAngles;
pos_ang = ang.*((ang>0)./(ang>0));
tt = 2*repmat(1:size(pos_ang,2),size(pos_ang,1),1);
scatter(tt(~isnan(pos_ang)),pos_ang(~isnan(pos_ang)),3,'k','filled','MarkerFaceAlpha',0.2,'MarkerEdgeAlpha',0.05)
hold on
ang = allJobs{trial_num}.Tracks.smoothAngles;
pos_ang = ang.*((ang>0)./(ang>0));
NN = 35;
clrs = distinguishable_colors(NN);
for kk = 1:NN
    plot(tt(1,:),pos_ang(kk,:),'Color',clrs(kk,:),'LineWidth',1.5)
end
set(gca,'XTickLabel','');
ylabel('\theta_{right}')

% For negative angles:
axes('Units','Centimeters','Position',[1.5 7.5 12 3])
ang = allJobs{trial_num}.Tracks.origAngles;
neg_ang = ang.*((ang<0)./(ang<0));
tt = 2*repmat(1:size(neg_ang,2),size(neg_ang,1),1);
scatter(tt(~isnan(neg_ang)),abs(neg_ang(~isnan(neg_ang))),3,'k','filled','MarkerFaceAlpha',0.2,'MarkerEdgeAlpha',0.05)
hold on
ang = allJobs{trial_num}.Tracks.smoothAngles;
neg_ang = ang.*((ang<0)./(ang<0));
clrs = distinguishable_colors(NN);
for kk = 1:NN
    plot(tt(1,:),abs(neg_ang(kk,:)),'Color',clrs(kk,:),'LineWidth',1.5)
end
ylabel('\theta_{left}')
xlabel('Time [ms]')



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

end
