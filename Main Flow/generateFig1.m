function [] = generateFig1(allJobsBeforeScaling,px2mm,vid_path,curv_example_job)
% Calls relevant subroutines:
[fig1_handle] = generateFig1VideoHeadExamples(allJobsBeforeScaling,px2mm,vid_path);
generateFig1AnglePanel(allJobsBeforeScaling,fig1_handle,vid_path,px2mm);
generateFig1CurvPanel(px2mm,fig1_handle,vid_path,curv_example_job);
generateFig1TraceExample(allJobsBeforeScaling,px2mm,fig1_handle);


% Fix positioning, 1cm down
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

end