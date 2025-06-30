function [legend_labels] = makePolarTicks(values)
% Turn the axis into nicely formatted polar values (instead of decimals)
[N,D] = rat(values/pi);
for i = 1:length(N)
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
    legend_labels{i} = txt_str;


end
end
