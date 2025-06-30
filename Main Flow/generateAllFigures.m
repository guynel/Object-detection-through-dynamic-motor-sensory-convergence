%{
% Add codes to path, replace with path to codes directory:
codes_path = ''
addpath(genpath(codes_path));

% Replace with path to vids directory
vid_path = '';

% Load data:
% Includes all the trials data, the scalar the converts pixels to mm, and
% data from a job that was eventually excluded from analyses, but used to
% visualize curvature:
% Replace with path:
data_path = '';
% Filename:
data_filename = 'Data.mat';
load(fullfile(data_path,data_filename));

allJobsBeforeScaling = allJobs;
%}
%% Plot Figure 1
generateFig1(allJobsBeforeScaling,px2mm,vid_path,curv_example_job);

%% Scale the distribution of angles in the theta-dot / kappa-dot plane
[allJobs,curv_scale_coeff] = scaleCurvature(allJobsBeforeScaling);

%% Generate figure 2
% Takes a few minutes to run if bootstrap tests are performed, see flag
% inside
% Also generates extended figure 1
[allJobs,cont_thresh_cm,alpha_star_val] = generateFig2(allJobs,px2mm,curv_scale_coeff);
cont_thresh_px = (cont_thresh_cm*10)*px2mm;
alpha_star_displacement = pi/2-alpha_star_val;

%% Adding cycle segmentation and approach
% Parse to cycles
[pump_struct,allJobs] = addCycleData(allJobs,cont_thresh_px);
% Define approach phase
extend_to_cycle_end_flag = false;
[allJobs] = defineApproachPhase(allJobs,px2mm,cont_thresh_px,extend_to_cycle_end_flag);
% In the cycle structure, add a marking of whether the cycle begins during the approach
[allJobs,pump_struct] = updatePumpStruct(allJobs,pump_struct);

%% Plot figure 3, extended figure 2
extend_to_cycle_end_flag = false;
[allJobs] = generateFig3(allJobs,px2mm,3,pi/2,extend_to_cycle_end_flag,cont_thresh_px,alpha_star_displacement);
generateSupp2(allJobs,px2mm,extend_to_cycle_end_flag,cont_thresh_px,alpha_star_val);

%% Generate figure 4, extended figure 3
extend_to_cycle_end_flag = false;
allJobs = generateFig4(allJobs,px2mm,extend_to_cycle_end_flag,cont_thresh_px,alpha_star_val);
generateSupp3(allJobs,px2mm,extend_to_cycle_end_flag,cont_thresh_px,alpha_star_val);

%%
extend_to_cycle_end_flag = true;
[allJobs] = generateFig5(allJobs,px2mm,pump_struct,extend_to_cycle_end_flag,alpha_star_displacement,cont_thresh_px);

%% Supps:
generateSupp5(allJobs);
generateSupp6(allJobs,px2mm);

