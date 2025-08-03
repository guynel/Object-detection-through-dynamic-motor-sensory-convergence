function [] = generateFig1(allJobsBeforeScaling,px2mm,vid_path,curv_example_job)
% Calls relevant subroutines:
[fig1_handle] = generateFig1VideoHeadExamples(allJobsBeforeScaling,px2mm,vid_path);
generateFig1AnglePanel(allJobsBeforeScaling,fig1_handle,vid_path,px2mm);
generateFig1CurvPanel(px2mm,fig1_handle,vid_path,curv_example_job);
generateFig1TraceExample(allJobsBeforeScaling,px2mm,fig1_handle);


%% Fix positioning, 1cm down
axesHandles = findall(gcf, 'Type', 'axes');
legendHandles = findall(gcf, 'Type', 'legend');
clrBars = findall(gcf, 'Type', 'Colorbar');
translate_down_x_cm = 1; 
for ax = axesHandles'
    ax.Units = 'centimeters';
    ax.Position(2) = ax.Position(2) - translate_down_x_cm;
end
for lgd = legendHandles'
    lgd.Units = 'centimeters';
    lgd.Position(2) = lgd.Position(2) - translate_down_x_cm;
end
for cbar = clrBars'
    cbar.Units = 'centimeters';
    cbar.Position(2) = cbar.Position(2) - translate_down_x_cm;
end
drawnow expose

%% Change dimensions of panels, to fit publication standards

% Get handles
allAxes = findall(gcf, 'type', 'axes');
allAxes = allAxes(~strcmp(get(allAxes, 'Tag'), 'Colorbar'));
for k = 1:length(allAxes)
    allAxes(k).Units = 'centimeters';
end
% Get the groups of axes by height-position
pos = cell2mat(get(allAxes, 'Position'));
pos = round(pos,2);
grp_idx = zeros(size(pos(:,2)));
grp_idx(pos(:,2) > 10) = 3;
grp_idx(pos(:,2) < 10) = 4;
grp_idx(pos(:,2) == 22.5) = 1;
grp_idx(pos(:,2) == 18.75) = 2;

%% Manually fix positions
group_one = find(grp_idx == 1);
[~,sorted_idx] = sort(pos(group_one,1));
group_one_idx = group_one(sorted_idx);
panel_size = 3.2;
for k = 1:length(group_one)
    allAxes(group_one_idx(k)).Position = [4+(k-1)*panel_size 23 panel_size panel_size]
end

group_two = find(grp_idx == 2);
[~,sorted_idx] = sort(pos(group_two,1));
group_two_idx = group_two(sorted_idx);
allAxes(group_two_idx(1)).Position = [3.8 20 2.5 2.5];
allAxes(group_two_idx(2)).Position = [7 20 2.5 2.5];
allAxes(group_two_idx(3)).Position = [8.64 20 2.5 2.5];
allAxes(group_two_idx(4)).Position = [10.5 20 6 2.5];
clr_bar = findall(gcf, 'type', 'Colorbar')
clr_bar.Position([2,4]) = allAxes(group_two_idx(4)).Position([2,4]);

leg = findall(gcf, 'type', 'legend');
leg.Position = [3.8 17.5 leg.Position(3:4)];
group_three = find(grp_idx == 3);
[~,sorted_idx] = sort(pos(group_three,1));
group_three_idx = group_three(sorted_idx);
allAxes(group_three_idx(1)).Position = [3.8 16.25 7.25 1.25];
set(allAxes(group_three_idx(1)),'XTickLabel',{})
allAxes(group_three_idx(2)).Position = [3.8 14.5 7.25 1.25];
allAxes(group_three_idx(3)).Position = [13.5 14.5 3 3];

group_four = find(grp_idx == 4);
[~,sorted_idx] = sort(pos(group_four,1));
group_four_idx = group_four(sorted_idx);
allAxes(group_four_idx(2)).Position = [3.8 11.25 7.25 1.25];
set(allAxes(group_four_idx(2)),'XTickLabel',{})
allAxes(group_four_idx(1)).Position = [3.8 9.5 7.25 1.25];
allAxes(group_four_idx(3)).Position = [13.5 11 1.5 1.5];
allAxes(group_four_idx(4)).Position = [13.5 9.5 1.5 1.5];
allAxes(group_four_idx(5)).Position = [15 11 1.5 1.5];
allAxes(group_four_idx(6)).Position = [15 9.5 1.5 1.5];
allAxes(group_four_idx(8)).Position = [16.3 11.2 1.15 1.15];
allAxes(group_four_idx(7)).Position = [16.3 9.6 1.15 1.15];
axes(allAxes(group_four_idx(7)));
axis(axis().*[0 1 0 1]+[-0.2 0 -0.2 0])
allTextHandles = findall(gcf, 'Type', 'text');
allTextHandles = allTextHandles([10,14,7,11]);
for k = 1:2
    allTextHandles(k).Position(1) = allTextHandles(k).Position(1)+0.6;
end
allTextHandles(1).String = {'Cont.';'whisker'};
allTextHandles(2).String = {'Non-cont.';'whisker'};
allTextHandles(3).String = {'Cont.'};
allTextHandles(4).String = {'Non-cont.'};

end