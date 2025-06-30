function [fig_handle] = makeFullPagePDF()
% Make a full-page sized figure, for plotting
fig_handle = figure;
pos = [0 0 21 29.7];
set(fig_handle, 'PaperUnits','centimeters');
set(fig_handle,'Units','Centimeters','Position',pos);
set(fig_handle, 'PaperSize', [pos(3) pos(4)]);
set(fig_handle,'PaperPosition',pos)

end