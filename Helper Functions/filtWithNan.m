function [output,Hd] = filtWithNan(signal,filt_data)
% Called by addCycleData, to generate lowpass and bandpass median angle
% signals; Rarely has missing data, so needs to be interpolated

output = NaN(size(signal));
for row_num = 1:size(signal,1) % Generalized, but actually typically just one row
    t_start = find(~isnan(signal(row_num,:)),1); % 1st non-nan
    t_end = find(~isnan(signal(row_num,:)),1,'last'); % last non-nan
    sig = signal(row_num,t_start:t_end); % Truncated signals

    containedNans = false;
    if any(isnan(sig)) % Interpolate if needed
        t = find(~isnan(sig));
        tt = find(isnan(sig));
        sig(tt) = interp1(t,sig(t),tt,'linear');
        containedNans = true;
    end
    try
        filtSignal = filtfilt(filt_data.sos, filt_data.g,sig); % Filter
        if containedNans
            filtSignal(tt) = NaN; % Revert to NaNs, don't want to rely on interpolations
        end
        output(row_num,t_start:t_end) = filtSignal;
    end
end

end