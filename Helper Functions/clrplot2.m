function [] = clrplot2(x,y,clrs,varargin)
% Used in generating Fig. 1 (panel c)
% plot x and y, with each coordinate colored by the color in clrs
if nargin >= 4
    LineStyle = varargin{1};
else
    LineStyle = '-';
end

% If we want to use scatter markers as well
scatterFlag = 0;
if length(LineStyle)==2
    if any(LineStyle ~= '--') & any(LineStyle ~= '-.')
        scatterFlag = 1;
    end
end
        
% Handle dimensions and add null third dimension for 'surface' function
x = x(:)';
y = y(:)';
z = zeros(size(x));

% Handle dimensions of clrs
orig_clrs = clrs;
clrs = permute(clrs,[3, 1, 2]);

% Re-arrange vectors
x = [x;x];
y = [y;y];
z = [z;z];
clrs = [clrs;clrs];

% If plotting also requires markers:
if scatterFlag
    if ismember(LineStyle(2),'osd<>^vph')
        scatter(x(1,:),y(1,:),[],orig_clrs,LineStyle(2),'filled');
    else
        scatter(x(1,:),y(1,:),[],orig_clrs,LineStyle(2));
    end
    hold on
    surface(x,y,z,clrs,...
        'facecol','no',...
        'edgecol','interp',...
        'LineStyle','-',...
        'linew',1.5);
% Just plotting
else
    surface(x,y,z,clrs,...
        'facecol','no',...
        'edgecol','interp',...
        'LineStyle',LineStyle,...
        'linew',1.5);
end
end
