% Add utils path
clc; clear; close all;  % Clear command window, workspace, and close all figures
addpath("./utils/");      % Add the utils directory to the search path
addpath("./group_analysis/");      % Add the directory to the search path
addpath("./fnc_computation/");      % Add the directory to the search path
%% Data Initialization
% Define the number of subjects, timepoints, components, and sampling time
num_subjects = 100;       % Number of subjects
num_timepoints = 150;   % Number of timepoints
num_components = 10;    % Number of components
TR = 2;                 % repetition time
bandwidth = [0.01 0.1];
% Generate example post-processed fMRI timecourses
subTcs = randn(num_subjects, num_timepoints, num_components);
subTcs = post_processing(subTcs, TR, bandwidth);

% Define groupings for subjects (example, 1 is schizophrenia, 2 is controls)
group = randi([1 2], num_subjects, 1);
age= randi([20 60], num_subjects, 1); % Generate random ages between 20-65
sex = randi([1 2], num_subjects, 1); % Male or female
site = randi([1 5], num_subjects, 1); 

%% Compute Correlation FNC
corr_fnc = compute_correlation_fnc(subTcs);

%% Compute DTW across subjects
win_size = calculate_window_size(TR, bandwidth);
[dtw_stand, dtw_norm] = compute_DTW(subTcs, win_size);

%% Convert DTW distances to FNC
dtw_stand_fnc = -zscore(dtw_stand, 0, "all");
dtw_norm_fnc = -zscore(dtw_norm, 0, "all");

%% Perform Wald test from GLM
covariates = table(age, sex, site); % Create covariates table
%Group analysis for correlation
[corr_q_values, corr_p_values] = compute_glm_group_difference(corr_fnc, group, covariates);

%Group analysis for standard DTW
[dtw_stand_q_values, dtw_stand_p_values] = compute_glm_group_difference(dtw_stand_fnc, group, covariates);

%Group analysis for normalized DTW
[dtw_norm_q_values, dtw_norm_p_values] = compute_glm_group_difference(dtw_norm_fnc, group, covariates);