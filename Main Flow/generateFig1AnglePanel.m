function [] = generateFig1AnglePanel(allJobs,fig1_handle,vid_path,px2mm)
trial_num = 95; % Example trial
frame_num = 223; % Example frame
vid_name = allJobs{trial_num}.metadata.videoBasename;


% Make figure
main_fig = fig1_handle;
axes();

% Get video img
v = VideoReader(fullfile(vid_path,[vid_name,'.avi']));
img = double(read(v,frame_num))./255;
img = rgb2gray(img).*255;
mask = zeros(size(img, 1), size(img, 2));
mask(840:end,1:600) = ones;
img(find(mask)) = NaN;
img = uint8(inpaint_nans(img,4));
imshow(img)

% Get whiskers
job = allJobs{trial_num};
frame_idx = job.frameIndices == frame_num;
wh = job.results{frame_idx}.whiskers.plotWhiskers;
label = job.results{frame_idx}.whiskers.label;

% Roate the image for nicer presentation
angle = 45;  % rotation angle in degrees
cx = 425; cy = 488;  % center of rotation
rotated_image = imrotate(img, angle, 'bilinear');
theta = -deg2rad(angle);
rotation_matrix = [cos(theta), -sin(theta); sin(theta), cos(theta)];
imshow(rotated_image);
hold on;

% Plot head-center line
plot(655*[1,1],ylim,'w--','LineWidth',1.5)

% Plot colored-patches for angles
angs = round(job.results{frame_idx}.whiskers.theta(1:4));
clrs = [1, 1, 1; 1.0000    0.7812    0.4975; 0.9, 0.61, 0.39];
clrs = [clrs;clrs(end,:)];
x_target = 655*[1,1,1];
y_target = [975,970,960];
radius = [30,42,55]*2;     % Radius of the arc
circ_theta = angs([1,2,4]);
for k = 3:-1:1
    theta_start = 90;
    theta_end = 90 - circ_theta(k);
    theta = linspace(deg2rad(theta_start), deg2rad(theta_end), 50);
    x_arc = x_target(k) + radius(k) * cos(theta);
    y_arc = y_target(k) - radius(k) * sin(theta);
    x_fill = [x_target(k), x_arc, x_target(k)];
    y_fill = [y_target(k), y_arc, y_target(k)];
    hold on;
    fill(x_fill, y_fill, clrs(k,:), 'FaceAlpha', 0.7, 'EdgeColor', 'None'); % Transparent red wedge
    
end

% Plot whiskers + continuation to line
coords = cell(0);
for i = [1,2,4]
    x = wh{i};
    coords{i} = (rotation_matrix * (x - [cx, cy])')' + [645 950];
    
    x_target = 655;
    ang_target = 180-(90-angs(i));
    theta_rad = deg2rad(ang_target);
    y_target = coords{i}(1,2) + (x_target - coords{i}(1,1)) * tan(theta_rad);
    lin_interp(i,:) = [x_target,y_target];
    coords{i} = [lin_interp(i,:);coords{i}];
    plot(coords{i}(:,1), coords{i}(:,2),':', 'Color',clrs(i,:), 'LineWidth', 1.5);
    
end

%axis([560 1200 588 1150])
% Place angle annotations
idx_shift = [200 200 200 150];
y_shift = [-40 -10 80 90]-30;
x_shift = [65 60 40 100]-100;
for k = [1,2,4]
    text(coords{k}(end-idx_shift(k),1)+x_shift(k),coords{k}(end-idx_shift(k),2)+y_shift(k),['\theta = ',num2str(angs(k)),char(176)],'Color',clrs(k,:),'FontSize',10,'FontWeight','normal');
    hold on
end


% For the whiskers we use a slightly different value, because they are at a
% lower plane; here we use the one that best aligns with the top of the
% object, to avoid introducing confusion. Eventually px2mm is mostly used
% for calculating curvature - which is re-scaled anyway when alpha is
% considerd; and contact - which is data-driven
one_cm_length_in_px = 10*px2mm;

plot(700+[0,one_cm_length_in_px],[750 750],'w','LineWidth',2);
text(700,700,'1cm','Color','w','FontSize',10,'FontWeight','normal')
set(gca,'Units','centimeters','Position',[1.5 19 3 3]);
axis([650 1100 650 1150])
set(gca,'Position',[2.4 19.75 2.85 2.85])

%%
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