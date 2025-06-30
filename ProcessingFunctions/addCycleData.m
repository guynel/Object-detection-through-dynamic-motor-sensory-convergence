function [pump_struct,allJobs] = addCycleData(allJobs,cont_thresh)
% Function used for cycle segmentation of trials

% Do we want to see how well it works?
visFlag = 0;

%% Design lowpass and high-pass filters

min_frames_cont = 3;
lp_freq = 80;
hp_freq = 7;
hp_filt = struct();
lp_filt = struct();
[z,p,k] = butter(6,hp_freq/500,'high');
[hp_filt.sos,hp_filt.g] = zp2sos(z,p,k); 
[z,p,k] = butter(6,lp_freq/500,'low');
[lp_filt.sos,lp_filt.g] = zp2sos(z,p,k); 


%% Main loop, iterate over trials

% Parameters:
min_trace_len = 0; % Too short to be considered a cycle
max_frames_interp = Inf; % How many interpolated frames disqualify a cycle
max_interp_prct = Inf; % What percent of interpolation disqualifies a cycle?
max_cycle_dist = 10; % If there is missing data, what's the furthest one cycle can start after another's end to call them consecutive?

for trial_num = 1:length(allJobs)
    for side_num = 1:2
        allJobs{trial_num}.Tracks.cycles_processed(side_num) = 0;

        % Some trials are not well-tracked enough
        not_nan_frames = sum(~isnan(allJobs{trial_num}.Tracks.angles))>1;
        if (sum(not_nan_frames)./length(not_nan_frames))<0.1
            continue
        end

        % Actual work:
        allJobs = getAnglePhase(allJobs,trial_num,side_num,hp_filt,lp_filt);
        allJobs = getCycleFromPhase(allJobs,trial_num,side_num);
        if visFlag
            visualizeCycles(allJobs,trial_num,side_num);
            pause
            clf
        end
        allJobs{trial_num}.Tracks.cycles_processed(side_num) = 1;
    end
end

% Finally:
% Get a different structure with all cycles
pump_struct = generatePumpStruct(allJobs);
% Categorize them 
pump_struct = markCycles(pump_struct, allJobs,cont_thresh,min_frames_cont);
% Also jot down categroies of preceding and succeeding cycles
pump_struct = relateConsecutiveCycles(pump_struct,min_trace_len,max_frames_interp,max_interp_prct,max_cycle_dist);
end

%%
function [allJobs] = getAnglePhase(allJobs,trial_num,side_num,hp_filt,lp_filt)
% Performing Hilbert transform to acquire signal phase

% Get raw angles
side = mode(allJobs{trial_num}.Tracks.side,2) == side_num;
side = side./side;
ang = abs(allJobs{trial_num}.Tracks.angles.*side);

% To aggregate over whiskers, use the integral of the median of the
% derivative
d_ang = multiStepDiff(ang,1);
med_d_ang = nanmedian(d_ang,1);
nan_idx = isnan(med_d_ang);
med_ang = cumsum(med_d_ang,'omitnan');
med_ang = med_ang-nanmean(med_ang(:));
med_ang = med_ang+nanmean(abs(ang(:)));

% Filter to get bandpass version
lp_med_ang = filtWithNan(med_ang,lp_filt);
bp_med_ang = filtWithNan(lp_med_ang,hp_filt);
nan_vec = sum(~isnan(ang),1)>1;
nan_vec = nan_vec./nan_vec;

% Get phase-signal
phase_sig = NaN(size(bp_med_ang));
if ~all(isnan(bp_med_ang))
    nan_idx = find(isnan(bp_med_ang));
    not_nan_idx = find(~isnan(bp_med_ang));
    bp_med_ang(nan_idx) = interp1(not_nan_idx,bp_med_ang(not_nan_idx),nan_idx);
    t = find(~isnan(bp_med_ang),1):find(~isnan(bp_med_ang),1,'last');
    phase_sig(t) =  angle(hilbert(bp_med_ang(t)));
end

% Assign results to return
allJobs{trial_num}.Tracks.cycles{side_num}.med_ang = med_ang;
allJobs{trial_num}.Tracks.cycles{side_num}.lp_med_ang = lp_med_ang;
allJobs{trial_num}.Tracks.cycles{side_num}.bp_med_ang = bp_med_ang;
allJobs{trial_num}.Tracks.cycles{side_num}.phase_sig = phase_sig;
allJobs{trial_num}.Tracks.cycles{side_num}.nan_vec = nan_vec;


end

function [allJobs] = getCycleFromPhase(allJobs,trial_num,side_num)

% Align phase-signal cycles to nearest peaks in the original data
phase_sig = allJobs{trial_num}.Tracks.cycles{side_num}.phase_sig;
lp_med_ang = allJobs{trial_num}.Tracks.cycles{side_num}.lp_med_ang;
[peak_vals,peak_idx] = findpeaks(phase_sig);
peak_idx(abs(peak_vals-pi)>(pi/8)) = [];
[trough_vals,trough_idx] = findpeaks(-1*phase_sig);
trough_idx(abs(trough_vals-pi)>(pi/8)) = [];
peak_idx = [peak_idx;ones(size(peak_idx))];
trough_idx = [trough_idx;-1*ones(size(trough_idx))];
peaks_and_troughs = [peak_idx, trough_idx];
[~,sort_idx] = sort(peaks_and_troughs(1,:));
peaks_and_troughs = peaks_and_troughs(:,sort_idx);
extremum_class = [peaks_and_troughs(2,1:end-1);peaks_and_troughs(2,2:end)]';
extremum_class = 1+find(all(extremum_class == [1,-1],2));
cycle_starts = peaks_and_troughs(1,extremum_class);
first_cycle_start = peaks_and_troughs(1,find(peaks_and_troughs(2,:) == -1,1));
if first_cycle_start ~= cycle_starts(1)
    cycle_starts = [first_cycle_start, cycle_starts];
end

peak_idx = -1+find((multiStepDiff(sign(multiStepDiff(lp_med_ang,1)),1))>0);
[~,idx] = min(abs(peak_idx'-cycle_starts));
if length(unique(idx))<length(cycle_starts)
    dup_idx = find(idx(2:end) == idx(1:end-1));
    idx(dup_idx) = [];
    cycle_starts(dup_idx) = [];
end
cycle_starts = peak_idx(idx);

% Arrange extermum in pairs to get cycle edges:
complete_cycles = ones(size(cycle_starts));
if cycle_starts(1) ~= 1
    cycle_starts = [1,cycle_starts];
    complete_cycles = [0, complete_cycles];
end
if cycle_starts(end) < length(phase_sig)
    cycle_starts = [cycle_starts, length(phase_sig)];
    complete_cycles = [complete_cycles, 0];
end
% Result
cycle_start_end = [cycle_starts(1:end-1); -1+cycle_starts(2:end)]';

% From the low-pass (not band-pass!) get the up and down traces of whisking
sign_vec = [];
for kk = [-1,1]
    kk_idx = double(kk == 1)+1;
    sign_vec(kk_idx,:) = sign(multiStepDiff(lp_med_ang,1)) == kk;
end
sign_vec = sign_vec./sign_vec;
sign_vec = [sign_vec(:,2:end), NaN(2,1)];

% How much of a specific cylce was intrepolated for the hilbert transform?
nan_vec = allJobs{trial_num}.Tracks.cycles{side_num}.nan_vec;
cycle_t = 1:size(nan_vec,2);
cycle_t = cycle_t>=cycle_start_end(:,1) & cycle_t<=cycle_start_end(:,2);
cycle_non_interp = (nan_vec.*cycle_t) == 1;
non_interp_prct = sum(cycle_non_interp,2)./sum(cycle_t,2);

% Assign results to return
cycles = allJobs{trial_num}.Tracks.cycles{side_num};
cycles.cycle_start_end = cycle_start_end;
cycles.complete_cycles = complete_cycles;
cycles.cycle_t = cycle_t;
cycles.cycle_non_interp = cycle_non_interp;
cycles.non_interp_prct = non_interp_prct;
cycles.sign_vec = sign_vec;
allJobs{trial_num}.Tracks.cycles{side_num} = cycles;

end

function [] = visualizeCycles(allJobs,trial_num,side_num)
% Visualizing for debugging

side = mode(allJobs{trial_num}.Tracks.side,2) == side_num;
side = side./side;
ang = abs(allJobs{trial_num}.Tracks.angles.*side);
nan_vec = allJobs{trial_num}.Tracks.cycles{side_num}.nan_vec;
med_ang = allJobs{trial_num}.Tracks.cycles{side_num}.med_ang;
lp_med_ang = allJobs{trial_num}.Tracks.cycles{side_num}.lp_med_ang;
bp_med_ang = allJobs{trial_num}.Tracks.cycles{side_num}.bp_med_ang;
phase_sig = allJobs{trial_num}.Tracks.cycles{side_num}.phase_sig;
sign_vec = allJobs{trial_num}.Tracks.cycles{side_num}.sign_vec;
cycle_start_end = allJobs{trial_num}.Tracks.cycles{side_num}.cycle_start_end;
complete_cycles = allJobs{trial_num}.Tracks.cycles{side_num}.complete_cycles;



subplot(4,1,1)
plot((ang.*nan_vec)','k');
hold on
plot(med_ang,'b','LineWidth',2)
plot(lp_med_ang,'r','LineWidth',2)
subplot(4,1,2);
hold on
plot(lp_med_ang,'k');
lp_med_ang = lp_med_ang.*nan_vec;
clrs = 'br';
for kk = 1:2
    plot(lp_med_ang.*abs(sign_vec(kk,:)),'-','Color',clrs(kk),'LineWidth',2);
    hold on
end
subplot(4,1,3);
plot(phase_sig);
ylim([-1 1]*pi);

for panel_num = 1:3
    subplot(4,1,panel_num)
    y_lim = ylim();
    y_lim(2) = y_lim(2)-y_lim(1);
    for cycle_num = 1:size(cycle_start_end,1)
        tt = cycle_start_end(cycle_num,:);
        tt(2) = tt(2)-tt(1);
        if complete_cycles(cycle_num)
            f_clr = [0 0 0];
            rectangle('Position',[tt(1),y_lim(1),tt(2),y_lim(2)],'FaceColor',[f_clr 0.1]);
        end        
    end
    hold on
end


end


function [pump_struct] = generatePumpStruct(allJobs)
% Extract data into a struct that allows easy work on cycles, without
% trials
counter = 1;
for trial_num = 1:length(allJobs)
    for side_num = 1:2
        if allJobs{trial_num}.Tracks.cycles_processed(side_num)
            for cycle_num = 1:size(allJobs{trial_num}.Tracks.cycles{side_num}.cycle_start_end,1)
                pump_struct(counter).trial_num = trial_num;
                pump_struct(counter).side_num = side_num;
                
                cycles = allJobs{trial_num}.Tracks.cycles{side_num};
                pump_struct(counter).complete_cycle = cycles.complete_cycles(cycle_num);
                pump_struct(counter).start_end = cycles.cycle_start_end(cycle_num,:);
                t = cycles.cycle_t(cycle_num,:);
                t = t./t;
                pump_struct(counter).t = t;
                pump_struct(counter).nan_vals = (cycles.nan_vec).*t;
                non_interp_vals = cycles.cycle_non_interp(cycle_num,:);
                non_interp_vals = (non_interp_vals./non_interp_vals).*t;
                pump_struct(counter).non_interp_vals = non_interp_vals;
                sign_vec = cycles.sign_vec;
                pump_struct(counter).prot_retract = sign_vec.*t;
                pump_struct(counter).non_interp_prct = cycles.non_interp_prct(cycle_num);

                pump_struct(counter).phase_sig = cycles.phase_sig.*t;
                pump_struct(counter).med_ang = cycles.med_ang.*t;
                pump_struct(counter).lp_med_ang = cycles.lp_med_ang.*t;
                pump_struct(counter).bp_med_ang = cycles.bp_med_ang.*t;
                pump_struct(counter).serial_num = counter;
                counter = counter+1;
                
            end
        end
    end
end
end


function [pump_struct] = markCycles(pump_struct,allJobs,cont_thresh,min_frames_cont)
% Labels by number of pumps, contact, gather individual whisker data
for pump_num = 1:length(pump_struct)
    trial_num = pump_struct(pump_num).trial_num;
    side_num = pump_struct(pump_num).side_num;
    
    % How many pumps
    islands = getOneIslands(pump_struct(pump_num).prot_retract == 1);
    islands = [[islands{1},-1*ones(size(islands{1},1),1)];...
        [islands{2},ones(size(islands{2},1),1)]];
    [~,trace_idx] = sort(islands(:,1));
    islands = islands(trace_idx,:);
    t = 1:size(pump_struct(pump_num).prot_retract,2);
    signed_traces = t >= islands(:,1) & t <= islands(:,2);
    pump_struct(pump_num).signed_traces = signed_traces./signed_traces;
    pump_struct(pump_num).trace_signs = islands(:,3);
    pump_struct(pump_num).prot_num = sum(pump_struct(pump_num).trace_signs == 1);
    pump_struct(pump_num).is_double_pump = pump_struct(pump_num).prot_num>1;
    
    % Is there contact
    t = pump_struct(pump_num).t;
    side = mode(allJobs{trial_num}.Tracks.side,2) == side_num;
    side = side./side;
    wh_obj_dist = allJobs{trial_num}.Tracks.whiskerObjDist;
    wh_obj_dist = wh_obj_dist.*side;
    cont = wh_obj_dist<cont_thresh;
    cont = cont./cont;
    pump_cont = cont.*t;
    pump_struct(pump_num).is_cont = sum(sum(pump_cont==1,1)) >= min_frames_cont;

    % Does it qualify as a tip
    if sum(pump_struct(pump_num).trace_signs == 1)<2
        pump_struct(pump_num).is_tip = 0;
    else
        prot_idx = find(pump_struct(pump_num).trace_signs == 1);
        prot_idx = prot_idx(end);
        last_prot_start = find(pump_struct(pump_num).signed_traces(prot_idx,:) == 1,1);
        tt = t;
        tt(last_prot_start:end) = NaN;
        contact_before_final_protraction = pump_cont.*tt;
        if any(contact_before_final_protraction(:) == 1)
            pump_struct(pump_num).is_tip = 1;
        else
            pump_struct(pump_num).is_tip = 0;
        end
    end

    ang = allJobs{trial_num}.Tracks.smoothAngles.*side.*t;
    curv = allJobs{trial_num}.Tracks.smoothCurv.*side.*t;
    d_ang = allJobs{trial_num}.Tracks.diffAngles.*side.*t;
    d_curv = allJobs{trial_num}.Tracks.diffCurv.*side.*t;
    nan_tracks = all(isnan(ang),2) & all(isnan(curv),2);
    pump_struct(pump_num).ang = ang(~nan_tracks,:);
    pump_struct(pump_num).curv = curv(~nan_tracks,:);
    pump_struct(pump_num).d_ang = d_ang(~nan_tracks,:);
    pump_struct(pump_num).d_curv = d_curv(~nan_tracks,:);
    pump_struct(pump_num).wh_obj_dist = wh_obj_dist(~nan_tracks,:);



    pump_struct(pump_num).snoutObjDist = allJobs{trial_num}.Tracks.snoutObjDist.*t;
    pump_struct(pump_num).homing_ang = allJobs{trial_num}.Tracks.homing_ang.*t;
    pump_struct(pump_num).head_ang = allJobs{trial_num}.Tracks.head_ang.*t;

end
end


function [pump_struct] = relateConsecutiveCycles(pump_struct,min_trace_len,max_frames_interp,max_interp_prct,max_cycle_dist)
% Denote categories not only for each cycle, but for its preceding and
% following cycles
rmv_idx = zeros(1,length(pump_struct));

%  Mark cycles for removal if:
% They have basically no data
% They are very short 
% They are over-interpolated
for pump_num = 1:length(pump_struct)
    if isempty(pump_struct(pump_num).ang)
        rmv_idx(pump_num) = 1;
    elseif sum(all(isnan(pump_struct(pump_num).ang),1)<min_trace_len)
        rmv_idx(pump_num) = 1;
    end

    non_interp_prct = pump_struct(pump_num).non_interp_prct;
    interp_prct = (1-non_interp_prct);
    cycle_len = 1+diff(pump_struct(pump_num).start_end);
    if cycle_len < min_trace_len
        rmv_idx(pump_num) = 1;
    end
    if interp_prct > max_interp_prct
        rmv_idx(pump_num) = 1;
    end
    if interp_prct*cycle_len > max_frames_interp
        rmv_idx(pump_num) = 1;
    end
end

% Actual removal
pump_struct(logical(rmv_idx)) = [];
for pump_num = 1:length(pump_struct)
    pump_struct(pump_num).serial_num = pump_num;
end

% For each cycle, find predecessor, if exists
for pump_num = 1:length(pump_struct)
    trial_num = pump_struct(pump_num).trial_num;
    side_num = pump_struct(pump_num).side_num;

    trial_side_idx = find([pump_struct.trial_num] == trial_num & ...
        [pump_struct.side_num] == side_num);
    trial_side_data = pump_struct(trial_side_idx);
    
    curr_cycle_start = pump_struct(pump_num).start_end(1);
    cycles_start_end = [trial_side_data.start_end];
    cycles_start_end = reshape(cycles_start_end,2,length(cycles_start_end )/2)';
    cycles_end = cycles_start_end(:,2);
    dist_from_cycle_start = cycles_end-curr_cycle_start;
    dist_from_cycle_start(dist_from_cycle_start>=0) = Inf;
    [min_dist_val,min_dist_idx] = min(abs(dist_from_cycle_start));
    if abs(min_dist_val) < max_cycle_dist
        pump_struct(pump_num).prev_cycle_serial = ...
            trial_side_data(min_dist_idx).serial_num; 
    else
        pump_struct(pump_num).prev_cycle_serial = NaN;
    end
end

% Find next cycle
for pump_num = 1:length(pump_struct)
    pump_struct(pump_num).next_cycle_serial = NaN;
end
for pump_num = 1:length(pump_struct)
    all_serial = [pump_struct.serial_num];
    prev_cycle_serial = pump_struct(pump_num).prev_cycle_serial;
    if ~isnan(prev_cycle_serial)
        idx = find(all_serial == prev_cycle_serial);
        pump_struct(idx).next_cycle_serial =  pump_struct(pump_num).serial_num;
    end  
end

% Gather conscutive cycle characteristics
for pump_num = 1:length(pump_struct)
    all_serial = [pump_struct.serial_num];
    next_cycle_serial = pump_struct(pump_num).next_cycle_serial;
    prev_cycle_serial = pump_struct(pump_num).prev_cycle_serial;

    if ~isnan(prev_cycle_serial)
        prev_cycle_idx = find(all_serial == prev_cycle_serial);
        pump_struct(pump_num).prev_is_cont = ...
            pump_struct(prev_cycle_idx).is_cont;
        pump_struct(pump_num).prev_is_double_pump = ...
            pump_struct(prev_cycle_idx).is_double_pump;
        pump_struct(pump_num).prev_is_tip = ...
            pump_struct(prev_cycle_idx).is_tip;
    else
        pump_struct(pump_num).prev_is_cont = NaN;
        pump_struct(pump_num).prev_is_double_pump = NaN;
        pump_struct(pump_num).prev_is_tip = NaN;

    end
    if ~isnan(next_cycle_serial)
        next_cycle_idx = find(all_serial == next_cycle_serial);
        pump_struct(pump_num).next_is_cont = ...
            pump_struct(next_cycle_idx).is_cont;
        pump_struct(pump_num).next_is_double_pump = ...
            pump_struct(next_cycle_idx).is_double_pump;
        pump_struct(pump_num).next_is_tip = ...
            pump_struct(next_cycle_idx).is_tip;
    else
        pump_struct(pump_num).next_is_cont = NaN;
        pump_struct(pump_num).next_is_double_pump = NaN;
        pump_struct(pump_num).next_is_tip = NaN;
    end

end
end