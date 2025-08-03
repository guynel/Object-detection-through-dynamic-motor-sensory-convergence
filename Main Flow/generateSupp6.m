function [] = generateSupp6(allJobs,px2mm)
%% Extended Fig. 6

% Show that head-position on the y-axis align nicely when plotted with
% respect to end of approach phase
makeFullPagePDF();
axes('Units','centimeters','Position',[1.5 15 12 3])
for trial_num = 1:length(allJobs)
    hc = 0.1*allJobs{trial_num}.Tracks.headCenter(2,:)./px2mm; % Pixels to cm
    trim_idx = allJobs{trial_num}.Tracks.trim_idx;
    if isnan(trim_idx)
        continue
    end
    tt = 2*((1:length(hc))-trim_idx); % Frames to ms
    plot(tt,-1*(hc-hc(trim_idx)),'Color',[0 0 0 0.1]);
    hold on
end
axis([-750 2000 -6 6 ]);
ylabel('Y_{head} [cm]');
xlabel('Time [ms]')



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

end