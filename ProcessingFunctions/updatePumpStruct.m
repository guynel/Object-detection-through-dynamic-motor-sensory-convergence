function [allJobs,pump_struct] = updatePumpStruct(allJobs,pump_struct)
% After approach phase has been defined, if a cycle happens during
% approach, label it as such
[pump_struct.during_approach] = deal(0);

for trial_num = 1:length(allJobs)
    trim_idx = allJobs{trial_num}.Tracks.trim_idx;
    if ~isnan(trim_idx)
        pump_idx = find([pump_struct.trial_num] == trial_num);
        for k = pump_idx
            if pump_struct(k).start_end(1)<trim_idx
                pump_struct(k).during_approach = 1;
            end
        end
    end

end


end