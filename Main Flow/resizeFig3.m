%%
% Get handles
allAxes = findall(gcf, 'type', 'axes');
allAxes = allAxes(~strcmp(get(allAxes, 'Tag'), 'Colorbar'));
allAxes = flipud(allAxes);
polarAxes = findall(gcf, 'type', 'polaraxes');
all_leg = findall(gcf, 'type', 'legend');

%% Manually fix positions, top half
allAxes(1).Position = [3 24.75 6 0.55];
ylbl = allAxes.YLabel;
allAxes(1).YColor = 'none';
allAxes(1).YLabel.Color = 'k';
allAxes(1).YLabel.Units = 'centimeters';
allAxes(1).YLabel.Position = [-1.3 0.2 1];
allAxes(1).YLabel.Rotation = 0;
allAxes(2).Position = [3 23.55 6 0.95];
allAxes(2).YLabel.Units = 'centimeters';
allAxes(2).YLabel.Position = [-1.3 0.75 1];
allAxes(2).YLabel.String =   {'\textsf{Mobility}';'$median(|\Delta{\alpha}|)$'};
allAxes(2).YLabel.Rotation = 0;
set(allAxes(2),'YAxisLocation','right')
allAxes(3).YLabel.Units = 'centimeters';
set(allAxes(3),'YAxisLocation','right')
allAxes(3).TickLabelInterpreter = 'latex';
allAxes(3).YLabel.Position = [-1.3 1 1];
allAxes(3).YLabel.String =   {'\textsf{Directional}';'\textsf{Tendency}';...
    '$p(\Delta\alpha > 0)-0.5$'};
allAxes(3).YLabel.Rotation = 0;
allAxes(3).Position = [3 22.1 6 0.95];
allAxes(4).Position = [11.25 21.75 4 4];
allAxes(5).Position = [3 22.1 6 3.2];
leg = all_leg(end);
leg.Position = [3 26.05 6 0.56];

%% Manually fix position, bottom half
h = findall(gcf,'Type','Histogram');
max_r = max([h.Values]);
max_r = 0.05+ceil(max_r*100)/100;
for k = 1:length(polarAxes)
    polarAxes(k).RLim = [0,max_r];
end
polarAxes = flipud(polarAxes);

polarAxes(1).Position = [1 17.85 1.85 1.85];
polarAxes(2).Position = [4.25 17.85 1.85 1.85];
polarAxes(3).Position = [7.25 17.85 1.85 1.85];
polarAxes(4).Position = [1 15.45 1.85 1.85];
polarAxes(5).Position = [4.25 15.45 1.85 1.85];
polarAxes(6).Position = [7.25 15.45 1.85 1.85];
polarAxes(7).Position = [1 13.05 1.85 1.85];
polarAxes(8).Position = [4.25 13.05 1.85 1.85];
polarAxes(9).Position = [7.25 13.05 1.85 1.85];
polarAxes(10).Position = [1 9.8 1.85 1.85];
polarAxes(11).Position = [4.25 9.8 1.85 1.85];
polarAxes(12).Position = [7.25 9.8 1.85 1.85];

% Now non-polar
allAxes(6).Position = [11 17.85 1.85 1.85];
allAxes(7).Position = [11 15.45 1.85 1.85];
allAxes(8).Position = [11 13.05 1.85 1.85];
allAxes(9).Position = [11 9.8 1.85 1.85];

for k = 1:(length(all_leg))
    all_leg(k).ItemTokenSize = [18,18];
end

all_leg(4).FontSize = 8;
all_leg(4).Position = [12.85 18.1 3 1];
all_leg(3).FontSize = 8;
all_leg(3).Position = [12.85 15.8 3 1];
all_leg(2).FontSize = 8;
all_leg(2).Position = [12.85 13.4 3 1];
all_leg(1).FontSize = 8;
all_leg(1).Position = [12.85 10.1 3 1];

%% Get text and re-position
allTextHandles = findall(gcf, 'Type', 'text');
allTextHandles = [allTextHandles;findall(gcf, 'Type', 'textbox')];
desired_text = {'First';'Interim';'Last';'Contact';'Same side';'Other side'};
text_idx = [];
for k = 1:length(allTextHandles)
    if ~isempty(allTextHandles(k).String)
        if any(cellfun(@(x) all(strcmpi(allTextHandles(k).String,x)==1),desired_text))
            text_idx(end+1) = k;
        end
    end
end
allTextHandles = allTextHandles(text_idx);

for k = 1:3
    allTextHandles(k).Position = [0.85 2.45 1]
end
for k = 4:6
    allTextHandles(k).Position = [0.9 13.5+2.4*(k-4) 1 1]
    allTextHandles(k).VerticalAlignment = 'middle'
end



%% Cleaning up text formatting

child_list = get(gcf,'children');
for k = 1:length(child_list)
    try
        set(child_list(k),'FontSize',9,'FontWeight','Normal');
    end
end
txt_list = findall(gcf,'Type','Text');
for k = 1:length(txt_list)
    set(txt_list(k),'FontSize',9);
end
axes_list = findall(gcf,'Type','Axes');
for k = 1:length(axes_list)
    set(axes_list(k),'FontSize',9,'FontWeight','Normal');
    set(get(axes_list(k),'XLabel'),'FontSize',9,'FontWeight','Normal');
    set(get(axes_list(k),'YLabel'),'FontSize',9,'FontWeight','Normal');
end
annot_list = findall(gcf,'Type','textboxshape');
for k = 1:length(annot_list)
    set(annot_list(k),'FontSize',9);
end
drawnow expose
