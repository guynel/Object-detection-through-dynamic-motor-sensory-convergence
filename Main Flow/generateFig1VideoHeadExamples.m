function [main_fig] = generateFig1VideoHeadExamples(allJobs,px2mm,vid_path)
% Example of head and whisker progress over trial
trial_num = 65; % Trial num
plot_frames = [30 110 155 208]; % Frames to show
vid_name = allJobs{trial_num}.metadata.videoBasename;
v = VideoReader(fullfile(vid_path,[vid_name,'.avi']));
main_fig = makeFullPagePDF();

for k = 1:length(plot_frames) % Iterate over the four frames:
    img = double(read(v,plot_frames(k)))./255;
    img = rgb2gray(img).*255;

    % Inpaint to remove a piece of hardware that was below the platform, so
    % as not visually mislead readers regarding elements on the platform
    % RATS HAVE NO ACCESS TO THIS AND CAN'T SEE IT!
    mask = zeros(size(img, 1), size(img, 2));
    mask(900:975,1:600) = ones;
    mask(976:end,1:600) = ones;
    img(find(mask)) = NaN;
    img = inpaint_nans(img,4);
    img = uint8(img);

    % Show frame
    fh(k) = axes('Units','centimeters','Position',[2.6+(k-1)*3.6 23.5 3.6 3.6]);
    imshow(img)
    % Plot whiskers
    frame_idx = find(allJobs{trial_num}.frameIndices == plot_frames(k));
    wh = allJobs{trial_num}.results{frame_idx}.whiskers.plotWhiskers;
    label = allJobs{trial_num}.results{frame_idx}.whiskers.label;
    wh = wh(label<50); %% Avoid small length tracks, which are associated with large track-label numbers; in other words - bad traces
    hold on
    cellfun(@(x) plot(x(:,1),x(:,2),'w'),wh)

    hold on
    plot([0,px2mm*10]+100,[1 1]*950,'w-','LineWidth',1.5)
    if k == 1
        text(100,800,'1cm','Color','w','FontSize',10,'FontWeight','Normal');
    end
end

%% Tidy up
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