function [res] = getDataPerCycle(allJobs,px2mm,side_pref_thresh,pre_cont_thresh_coeff,cont_thresh_px,extend_to_cycle_end_flag)
% Used in Fig. 3-5; cycle-based data for easier processing

pre_cont_thresh = pre_cont_thresh_coeff*cont_thresh_px; % How far a whisker has to be to be called pre-contact; isn't really used

res = struct();

% Pre contact data
res.pre_contact = cell(0);
for trial_num = 1:length(allJobs) % Iterate over trials
 
    trim_idx = allJobs{trial_num}.Tracks.trim_idx; % End of approach
    % Verify validity of approach data
    if isnan(trim_idx)
        continue
    end
    if ~isfield(allJobs{trial_num}.Tracks,'cycles')
        continue
    end
    pre_contact = find(any(allJobs{trial_num}.Tracks.whiskerObjDist<pre_cont_thresh,1),1);
    % Get data for pre-contact
    % Pre-contact we can theoretically extract from any trial
    for side_num = 1:2
        side = mode(allJobs{trial_num}.Tracks.side,2) == side_num;
        side = side./side;
        cycle_t = allJobs{trial_num}.Tracks.cycles{side_num}.cycle_start_end;
        cycle_t = cycle_t(cycle_t(:,2)<=pre_contact,:);

        cycle_seq = cell(0);
        for cycle_num = 1:size(cycle_t,1)
            tt = cycle_t(cycle_num,:);
            tt = tt(1):tt(end);
            d_ang = side.*allJobs{trial_num}.Tracks.diffAngles(:,tt);
            d_ang = side.*d_ang;
            d_curv = allJobs{trial_num}.Tracks.diffCurv(:,tt);
            cycle_alpha = cart2pol(d_ang,d_curv);
            cycle_alpha = cycle_alpha(:);
            cycle_alpha = cycle_alpha(~isnan(cycle_alpha));
            cycle_seq{cycle_num} = cycle_alpha;
        end
        if ~isempty(cycle_seq)
            res.pre_contact{end+1} = cycle_seq;
        end
    end
end


% Contact and non-contact data
% Pre-allocate:
res.same_side = cell(0);
res.contact = cell(0);
res.r.contact = cell(0);
res.RD.contact = cell(0);
res.tt.contact = cell(0);
res.theta.contact = cell(0);
res.kappa.contact = cell(0);
res.other_side = cell(0);
res.shape = [];
res.texture = [];
res.trial_num = [];
res.side_num = [];

for trial_num = 1:length(allJobs) % Iterate over trials

    trim_idx = allJobs{trial_num}.Tracks.trim_idx; % End of approach
    if isnan(trim_idx)
        continue
    end
    
    % Find first contact:
    dd = allJobs{trial_num}.Tracks.whiskerObjDist(:,1:trim_idx)<cont_thresh_px;
    cont = dd./dd;
    [~,first_cont] = find(dd == 1,1);
    if isempty(first_cont)
        continue
    end
    non_cont = (~dd)./(~dd); % Non-contact data 
    % Contact and non-contact sides:
    side = mode(allJobs{trial_num}.Tracks.side,2);
    cont_side = side.*cont;
    % Only use trials which have pronounced preferred unilateral contact on
    % one side (defined by side_pref_thresh)
    cont_sides = histcounts(cont_side,0.5:2.5);
    [~,pref_side] = max(cont_sides);
    if (cont_sides(pref_side)./sum(cont_sides))<side_pref_thresh
        continue
    end

    % Define again, because other side may have cycles going over the trim
    % idx:
    dd = allJobs{trial_num}.Tracks.whiskerObjDist<cont_thresh_px;
    cont = dd./dd;
    non_cont = (~dd)./(~dd);

    RD = allJobs{trial_num}.Tracks.onWhiskerDist.rad_dist;
    RD = 0.1*RD/px2mm;

    % Cont side:
    side_num = pref_side;
    side = mode(allJobs{trial_num}.Tracks.side,2) == side_num;
    side = side./side;
    cycle_t = allJobs{trial_num}.Tracks.cycles{side_num}.cycle_start_end;
    cycle_t = cycle_t((cycle_t(:,1)<trim_idx) & (cycle_t(:,2)>first_cont),:);
    if ~extend_to_cycle_end_flag
        cycle_t(end) = trim_idx;
    end
    % Pre-allocate
    cycle_seq_cont = cell(0);
    cycle_seq_non_cont = cell(0);
    cycle_seq_cont_r = cell(0);
    cycle_seq_cont_RD = cell(0);
    cycle_seq_cont_theta = cell(0);
    cycle_seq_cont_kappa = cell(0);
    % Get data for contact-side
    for cycle_num = 1:size(cycle_t,1)
        tt = cycle_t(cycle_num,:);
        tt = tt(1):tt(end);
        d_ang = allJobs{trial_num}.Tracks.diffAngles(:,tt);
        d_ang = d_ang.*side.*cont(:,tt);
        d_curv = allJobs{trial_num}.Tracks.diffCurv(:,tt);
        [cycle_alpha,cycle_r] = cart2pol(d_ang,d_curv);
        cyc_d_ang = d_ang(~isnan(cycle_alpha));
        cyc_d_curv = d_curv(~isnan(cycle_alpha));
        cycle_alpha = cycle_alpha(:);
        cycle_alpha = cycle_alpha(~isnan(cycle_alpha));
        cycle_alpha_cont = cycle_alpha;
        cycle_r = cycle_r(:);
        cycle_r = cycle_r(~isnan(cycle_r));
        cycle_r_cont = cycle_r;

        
        cycle_RD = RD(:,tt).*side.*cont(:,tt);
        cycle_RD(isnan(d_ang)) = NaN;
        cycle_RD = cycle_RD(:);
        cycle_RD = cycle_RD(~isnan(cycle_RD));
        cycle_RD_cont = cycle_RD;

        d_ang = allJobs{trial_num}.Tracks.diffAngles(:,tt);
        d_ang = d_ang.*side.*non_cont(:,tt);
        d_curv = allJobs{trial_num}.Tracks.diffCurv(:,tt);
        cycle_alpha = cart2pol(d_ang,d_curv);
        cycle_alpha = cycle_alpha(:);
        cycle_alpha = cycle_alpha(~isnan(cycle_alpha));
        cycle_alpha_non_cont = cycle_alpha;

        cycle_seq_cont{cycle_num} = cycle_alpha_cont;
        cycle_seq_cont_r{cycle_num} = cycle_r_cont;
        cycle_seq_cont_RD{cycle_num} = cycle_RD_cont;
        cycle_seq_cont_theta{cycle_num} = cyc_d_ang;
        cycle_seq_cont_kappa{cycle_num} = cyc_d_curv;

        cycle_seq_non_cont{cycle_num} = cycle_alpha_non_cont;
        cycle_seq_non_cont{cycle_num} = cycle_alpha_non_cont;
    end
    % Assign if cycle is not empty of contact data
    if ~isempty(cycle_seq_cont)
        res.contact{end+1} = cycle_seq_cont;
        res.r.contact{end+1} = cycle_seq_cont_r;
        res.RD.contact{end+1} = cycle_seq_cont_RD;
        res.theta.contact{end+1} = cycle_seq_cont_theta;
        res.kappa.contact{end+1} = cycle_seq_cont_kappa;
        res.same_side{end+1} = cycle_seq_non_cont;
        res.tt.contact{end+1} = cycle_t;
        res.shape(end+1) = allJobs{trial_num}.metadata.shape;
        res.texture(end+1) = allJobs{trial_num}.metadata.texture;
        res.trial_num(end+1) = trial_num;
        res.side_num(end+1) = side_num;
    else
        continue
    end
    


    % Same: Non-cont side, just for alpha
    side_num = ~(pref_side-1)+1;
    side = mode(allJobs{trial_num}.Tracks.side,2) == side_num;
    side = side./side;
    cycle_t = allJobs{trial_num}.Tracks.cycles{side_num}.cycle_start_end;
    %cycle_t = allJobs{trial_num}.Tracks.cycles{~(side_num-1)+1}.cycle_start_end;
    cycle_t = cycle_t((cycle_t(:,1)<trim_idx) & (cycle_t(:,2)>first_cont),:);
    if ~extend_to_cycle_end_flag
        cycle_t(end) = trim_idx;
    end
    cycle_seq_other_side = cell(0);
    for cycle_num = 1:size(cycle_t,1)
        tt = cycle_t(cycle_num,:);
        tt = tt(1):tt(end);
        d_ang = allJobs{trial_num}.Tracks.diffAngles(:,tt);
        d_ang = d_ang.*side.*non_cont(:,tt);
        d_curv = allJobs{trial_num}.Tracks.diffCurv(:,tt);
        cycle_alpha = cart2pol(d_ang,d_curv);
        cycle_alpha = cycle_alpha(:);
        cycle_alpha = cycle_alpha(~isnan(cycle_alpha));
        cycle_alpha_other_side = cycle_alpha;

        cycle_seq_other_side{cycle_num} = cycle_alpha;
    end
    if ~isempty(cycle_seq_other_side)
        res.other_side{end+1} = cycle_seq_other_side;
    end
end

% Lose cycles before first contact on the contact side
for trial_num = 1:length(res.contact)
    cycle_seq = res.contact{trial_num};
    non_empty_cyc = cellfun(@(x) ~isempty(x),cycle_seq);
    keep_idx = find(non_empty_cyc,1):find(non_empty_cyc,1,'last');
    res.contact{trial_num} = cycle_seq(keep_idx);
    cycle_seq = res.r.contact{trial_num};
    res.r.contact{trial_num} = cycle_seq(keep_idx);
    cycle_seq = res.RD.contact{trial_num};
    res.RD.contact{trial_num} = cycle_seq(keep_idx);
    cycle_seq = res.tt.contact{trial_num};
    res.tt.contact{trial_num} = cycle_seq(keep_idx,:);
    cycle_seq = res.theta.contact{trial_num};
    res.theta.contact{trial_num} = cycle_seq(keep_idx);
    cycle_seq = res.kappa.contact{trial_num};
    res.kappa.contact{trial_num} = cycle_seq(keep_idx);

    cycle_seq = res.same_side{trial_num};
    res.same_side{trial_num} = cycle_seq(keep_idx);
end

% Rarely, we will have a non-contact cylce between contact-cycles; as its
% impossible to number contacts in this scenario, remove these
rmv_idx = zeros(1,length(res.contact));
for trial_num = 1:length(res.contact)
    cycle_seq = res.contact{trial_num};
    if length(cycle_seq) == 0
        rmv_idx(trial_num) = 1;
    end
    if any(cellfun(@(x) length(x), cycle_seq) == 0)
        rmv_idx(trial_num) = 1;
    end
end
rmv_idx = logical(rmv_idx);

res.contact = res.contact(~rmv_idx);
res.same_side = res.same_side(~rmv_idx);
res.r.contact = res.r.contact(~rmv_idx);
res.RD.contact = res.RD.contact(~rmv_idx);
res.tt.contact = res.tt.contact(~rmv_idx);
res.theta.contact = res.theta.contact(~rmv_idx);
res.kappa.contact = res.kappa.contact(~rmv_idx);
res.shape(rmv_idx) = [];
res.texture(rmv_idx) = [];
res.trial_num(rmv_idx) = [];
res.side_num(rmv_idx) = [];
res.other_side = res.other_side(~rmv_idx);



end