function [] = mockPolarAxes(numThetaWedges,rTicks,axesHandle,rScale,rTicksClr,varargin)
% For advanced plotting of polar data, where MATLAB's built-in doesn't cut
% it. User specific number of circle cut into wedges; you can then plot as
% you wish on it - for example patches - because axes is still technically
% Cartesian
if length(varargin) == 0
    exclude_r_ticks = 0;
else
    exclude_r_ticks = varargin{1};
end

set(gcf,'CurrentAxes',axesHandle);
hold on

thetaVals = linspace(-pi,pi,361);
[x,y] = pol2cart(thetaVals',rTicks);
plot(x,y,'Color',[1 1 1]*0.75)

thetaVals =  linspace(-pi,pi,numThetaWedges+1);
[x,y] = pol2cart(thetaVals,rTicks');
plot(x,y,'Color',[1 1 1]*0.75);
axis equal

set(gca,'Color','None','XColor','None','YColor','None')

if ~exclude_r_ticks
for i = 1:length(rTicks)
    [x,y] = pol2cart(deg2rad(80),rTicks(i));
    if length(varargin)>1
        r_txt = varargin{2}{i};
    else
        r_txt = num2str(rTicks(i));
    end
    if nargin == 4
        text(x,y,r_txt,'Color',[1 1 1]);
    else
        text(x,y,r_txt,'Color',rTicksClr);
    end
end
end

[N,D] = rat(thetaVals/pi);
for i = 2:length(N)
    [x,y] = pol2cart(thetaVals(i),rScale*rTicks(end));
    num = num2str(N(i));
    denom = num2str(D(i));
    
    if D(i) == 1
        if abs(N(i)) ~= 1
            txt_str = ['$',num,'\pi$'];
        elseif N(i) == 1
            txt_str = ['$\pi$'];
        elseif N(i) == -1
            txt_str = ['$-\pi$'];
        end
    else
        if N(i)<0
            txt_str = ['$-\frac{',num2str(abs(N(i))),'}{',denom,'}\pi$'];
        else
            txt_str = ['$\frac{',num2str(N(i)),'}{',denom,'}\pi$'];
        end
    end
    if N(i) == 0
        txt_str = '0';
    end
   
    hAlign = 'center';
    vAlign = 'middle';

    tt = text(x,y,txt_str,...
        'FontWeight','bold','FontSize',16,...
        'Interpreter','latex','VerticalAlignment',...
        vAlign,'HorizontalAlignment',hAlign);
    
end