function [] = generateFig1CurvPanel(px2mm,fig1_handle,vid_path,job)
% A nice example of curvature 

% Get example job; From a video that was eventually
% excluded from the data set, hence loaded here to serve only as an
% illustration

% Get example video
vid = 'curv_example_vid.avi';
v = VideoReader(fullfile(vid_path,vid));

px2m = 1000*px2mm; % Units of m^-1

%% Get maximal curvature across all whiskers, to find range
whLen = 0;
curv = nan(size(job.Tracks.angles));
for mm = 1:size(job.Tracks.angles,1)
    
    single_wh_curv = getCurv(job, mm, whLen, 'max',px2mm);
    curv(mm,:) = single_wh_curv;
end

% Some heuristics for a nice illustration
[~,nan_idx] = getOneIslands(isnan(curv),3,0,0);
curv(~~nan_idx) = NaN;
curv(find(abs(multiStepDiff(curv,1))>50)) = NaN; % Remove exterme jumps, for lowpass in the following line
curv = lowpassFilt(curv,'Fst',200); % Lowpass

%% 
% Make figure
figure(fig1_handle);
all_vals = [];
maxLen = 0;
for i = 1:length(job.results)
    for j = [7,10,40] % These are 3 tracks all describing the same whisker; this was manually corrected to show curvature over long time-scale
        label = job.results{i}.whiskers.label;
        if any(label == j)
            whCurv = job.results{i}.whiskers.curvAlongWhisker{label == j};
            d = job.results{i}.whiskers.plotWhiskers{label==j};
            d = cumsum([0;sqrt(sum(diff(d,1).^2,2))]);
            maxLen = max(maxLen,max(d));
            all_vals = [all_vals,whCurv];
        end
    end
end
all_vals = -px2m*all_vals; % When re-extracting curvature, the sign needs to be flipped to fit convention
% Attempt to find decent range of values
zs = (all_vals-mean(all_vals))/std(all_vals);
all_vals(abs(zs)>2) = []; 
vals = [min(all_vals),max(all_vals)];
vals = [-30 30]; % Eventually hand-picked to be round and nice
clrs = jet(round(diff(vals))+2); % Colorscale for curvature

% Collect all relevant whisker-trace data
all_wh = [];
for i = 1:length(job.results)
    
    
    for j = [7,10,40]
        label = job.results{i}.whiskers.label;
        if any(label == j)
                label = find(label == j);
                wh = job.results{i}.whiskers.plotWhiskers{label};
                wh = transAndRot(job.results{i}.alignmentMat.alignMat,wh);
                wh = wh-wh(1,:);
                wh(:,2) = wh(:,2);
                all_wh = [all_wh;wh];
        end
    end
    if i == length(job.results)
        lastWh = wh;
    end
end
spaceFactor = 2;
xMax = spaceFactor*length(job.results)+max(lastWh(:,1));
yMin = min(all_wh(:,2));
yMax = max(all_wh(:,2));
ax = [0 xMax yMin yMax];
scale = (ax(2)-ax(1))/(ax(4)-ax(3));

% Set up panel
h = gobjects(1,4);
h(1) = axes();
hold on
axis(ax);
h(1).Units = 'centimeters';
widthCm = 16;
h(1).InnerPosition = [1 6 widthCm widthCm/scale];
plot(NaN,NaN,'Parent',h(1));
set(h(1),'XTick',[],'YTick',[],'XTickLabel',[],'YTickLabel',[],'Box','off','xcolor','None','ycolor','None');
set(gca,'Color','w')

% Go over observations and plot
for i = 1:length(job.results)
    
    for j = [7,10,40]
        label = job.results{i}.whiskers.label;
        if any(label == j)
                label = find(label == j);
                wh = job.results{i}.whiskers.plotWhiskers{label};
                wh = transAndRot(job.results{i}.alignmentMat.alignMat,wh); % Rotate to head-aligned
                wh = wh-wh(1,:);
                wh(:,2) = wh(:,2);
                wh(:,1) = wh(:,1)+(i-1)*spaceFactor;

                kappa = job.results{i}.whiskers.curvAlongWhisker{label}; % Curvature along whisker
                kappa = -px2m*kappa; % Correct kappa value
                kappa(kappa>max(vals)) = max(vals); % Color-code is bounded to avoid extremes
                kappa(kappa<min(vals)) = min(vals); % Color-code is bounded to avoid extremes
                kappa = (kappa-min(vals))+1;
                kappa = clrs(round(kappa),:); % Get clrs

                clrplot2(wh(:,1),-wh(:,2),kappa); % Plot
                hold on
                ylim([-300 520])
        end
    end
end
plot(2*[346 346],ylim,'k') % Vertical line for reference video frame
axes(h(1));
% Tidy up panel
h(1).Position = [6.95 8.275 9 3.63];
plot(xlim,[0 0],'k--');
set(gca,'XColor','k','XTick',0:200:(600*2),'XTickLabel',0:100:(500*2));
set(gca,'XTickLabelRotation',0);
xlabel('Time [ms]');
colormap(gca,jet)
clr_bar = colorbar('east');
clr_bar.Units = 'centimeters';
clr_bar.Position = [19 19 0.3 3.6];
clr_bar.Label.String = '\kappa [m^{-1}]';
clr_bar.Label.Rotation = 0;
clr_bar.Limits = [0,1];
set(clr_bar,'YTick',[0 1],'YTickLabel',{num2str(min(vals));num2str(max(vals))});
set(clr_bar,'Color','k');
set(clr_bar,'TickLabelInterpreter','Latex');
clr_bar.Label.Rotation = 90;

% Add reference frame and plot all whiskers with colored curvature on it
N = 1+ceil(diff(vals)/2);
h(2) = axes();
rotation_angle = deg2rad(-1); % Slight additional correction to align traces, due to rounding errors
R = [cos(rotation_angle), -sin(rotation_angle);
sin(rotation_angle),  cos(rotation_angle)];
for i = 346 % Relevant frame number
    idx = find(cellfun(@(x) x.frameNum,job.results)==i); % Get frame idx in the result field
    imshow(imrotate(read(v,idx),30,'bilinear'),'Parent',h(2)); % rotate image
    hold on;
    
    % Whisker data
    wh = job.results{idx}.whiskers.plotWhiskers;
    
    for j = 1:length(wh) % Iterate over whiskers


        wh = job.results{i}.whiskers.plotWhiskers{j};
        if pdist2(wh(1,:),wh(end,:))<30 % Bad traces (this job was not fully processed, as eventually not included in dataset)
            continue
        end
        % Transform whisker traces to image coordinates
        wh = transAndRot(job.results{i}.alignmentMat.alignMat,wh); 
        wh = (R*wh')';
        wh = wh+[678,895];

        % Associate curvature corrected and bound as above to obtain a
        % color-code
        kappa = job.results{i}.whiskers.curvAlongWhisker{j};
        kappa = -px2m*kappa;
        kappa(kappa>max(vals)) = max(vals);
        kappa(kappa<min(vals)) = min(vals);
        kappa = (kappa-min(vals))+1;
        kappa = clrs(round(kappa),:);
        hold on
        clrplot2(wh(:,1),wh(:,2),kappa); % plot
        drawnow
    end
end
% Tidy frame
axis([630 1175 99 1169])
h(2).Units = 'centimeters';
h(2).Position = [1.8 6 4 5.9];
plot(xlim,927*[1,1],'w--');
one_cm_length_in_px = px2mm*10;
plot([700 700]+[0,one_cm_length_in_px],1080*[1,1],'w','LineWidth',1.5)
text(700,1120,'1cm','Color','w','FontSize',10,'FontWeight','Normal')
axis([630 1175 510 1169 -1 1])
h(2).Position = [2.05 8.25 3.6 3.6];

% Repeat the proces but only for the whisker of interest
h(3) = axes();
rotation_angle = deg2rad(-1);
R = [cos(rotation_angle), -sin(rotation_angle);
sin(rotation_angle),  cos(rotation_angle)];
for i = 346
    idx = find(cellfun(@(x) x.frameNum,job.results)==i);
    imshow(imrotate(read(v,idx),30,'bilinear'),'Parent',h(3));
    hold on;

    wh = job.results{idx}.whiskers.plotWhiskers;
    
    for j = 4

        wh = job.results{i}.whiskers.plotWhiskers{j};
        if pdist2(wh(1,:),wh(end,:))<30
            continue
        end
        wh = transAndRot(job.results{i}.alignmentMat.alignMat,wh);
        wh = (R*wh')';
        wh = wh+[678,895];

        kappa = job.results{i}.whiskers.curvAlongWhisker{j};
        kappa = -px2m*kappa;
        kappa(kappa>max(vals)) = max(vals);
        kappa(kappa<min(vals)) = min(vals);
        kappa = (kappa-min(vals))+1;
        kappa = clrs(round(kappa),:);
        hold on
        clrplot2(wh(:,1),wh(:,2),kappa);
        drawnow
    end
end
% Plot a nice osculating circle for intuition
max_idx = 685;
kappa = job.results{i}.whiskers.curvAlongWhisker{j};
dp = diff(wh(max_idx+[-1,1],:));
wh_tang = dp / norm(dp);
wh_normal = -1*sign(kappa(max_idx)).*([0 -1; 1 0]*wh_tang')';
circ_center = wh(max_idx,:)+(1/kappa(max_idx))*wh_normal;

theta = linspace(0, 2*pi, 100);
circ_coords = circ_center + (1/kappa(max_idx)) * [cos(theta') sin(theta')];
circ_handle = plot(circ_coords(:,1), circ_coords(:,2), 'w');
uistack(circ_handle,'down');
plot([circ_center(1),wh(max_idx,1)], [circ_center(1),wh(max_idx,2)], 'w');
scatter(circ_center(1),circ_center(1),'wo','filled');


% Tidy up panel
axis([700 1000 99 1169])
h(3).Units = 'centimeters';
h(3).Position = [5.1 6 2 5.9];
axes(h(2));
ax = axis;
axes(h(3));
axis([695 1000 ax(3:4)]);
h(3).Position = [4.3 8.25 3.6 3.6];
plot(xlim,927*[1,1],'w--');
text(770,820,'r','Color','w','FontSize',10);

% Correct panel sizes
basic_panel_size = 2.85;
h(1).Position = [9.5 19.75 7 basic_panel_size];
h(2).Position = [5.5 19.75 basic_panel_size basic_panel_size];
h(3).Position = [7.35 19.75 basic_panel_size basic_panel_size];
clr_bar.Position = [16.65 19.75 0.3 basic_panel_size];
clr_bar.Label.Position = clr_bar.Label.Position+[-1 0 0 ];
%% Tidy everything up

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