function [curv] = getCurv(job, whNum, whLen, curvType,px2mm)
% This goes over all whiskers observations associated with a specific
% whisker label in a given job, and returns one of a number of curvature
% parameters based on input
curv = NaN(1,size(job.Tracks.curv,2)); % Pre-allocated for output

for frameIdx = 1:length(job.results) % Go over the results field
    frameNum = job.results{frameIdx}.frameNum; % Video frame numbers
    wh_idx = find(job.results{frameIdx}.whiskers.label == whNum); % Find label-associated data
    if isempty(wh_idx)
        continue
    end
    whCurv = job.results{frameIdx}.whiskers.curvAlongWhisker{wh_idx};
    whCurv = -px2mm*1000*whCurv; % Convert to meter^-1; flip sign, because original convention was arbitrarily flipped
    if all(isnan(whCurv(:)))
        continue
    end
    whDist = job.results{frameIdx}.whiskers.plotWhiskers{wh_idx}; % Length on whiskers
    whDist = cumsum(sqrt(sum((diff(whDist,[],1)).^2,2))); % Cumulative
    [~,wh_len_idx] = min(abs(whDist-whLen)); % Trim to some length, if desired by user
    whCurv = whCurv(1:wh_len_idx); 
    switch curvType
        case 'max'
            curv(frameNum) = max(whCurv);
            %curv(frameNum) = prctile(whCurv,50);
        case 'min'
            curv(frameNum) = min(whCurv);
        case 'median'
            curv(frameNum) = nanmedian(whCurv);
    end
end

end