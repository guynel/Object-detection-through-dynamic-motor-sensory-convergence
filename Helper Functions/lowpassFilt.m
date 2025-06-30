function [output,Hd] = lowpassFilt(signal,varargin)
% Used in Fig. 4 panel c, and in head-variable smoothing in Fig. 4, also
% extended Fig. 2 and some behind the scene processing

% Lowpass filter data, can receive filter specification
output = NaN(size(signal));

fieldNames = {'Fp';'Fst';'Ap';'Ast';'fs'}; % Default is specification not given
fieldValues = {1;250;1;20;500};
% Generate filter opts
if length(varargin)>0
    inputs = struct(varargin{:});
    inputNames = fieldnames(inputs);
    for i = 1:length(inputNames)
        idx = find(strcmpi(inputNames{i},fieldNames));
        if ~isempty(idx)
            fieldValues{idx} = inputs.(inputNames{i});
        end
    end
end
temp = [fieldNames,fieldValues]';
opts = struct(temp{:});

% Generate filter
Hd = fdesign.lowpass('Fp,Fst,Ap,Ast',opts.Fp/opts.fs,opts.Fst/opts.fs,opts.Ap,opts.Ast);
Hd = design(Hd,'butter');
% For each row in the input, filter
for i = 1:size(signal,1)
    % Truncate singal to actual start and finish (without missing data)
    t_start = find(~isnan(signal(i,:)),1);
    t_end = find(~isnan(signal(i,:)),1,'last');
    sig = signal(i,t_start:t_end);

    % Interpolate if it has NaNs
    containedNans = false;
    if any(isnan(sig))
        t = find(~isnan(sig));
        tt = find(isnan(sig));
        sig(tt) = interp1(t,sig(t),tt,'linear');
        containedNans = true;
    end
    warning('off') % Annoying filtfilt issue
    try
        % Filter
        filtSignal = filtfilt(Hd.sosMatrix, Hd.scaleValues,sig);
        
        if containedNans
            filtSignal(tt) = NaN; % Get rid of interpolated data
        end
        output(i,t_start:t_end) = filtSignal; % Output
    end
    warning('on');
end

end