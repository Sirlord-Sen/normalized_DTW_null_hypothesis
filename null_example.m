% Add utils path
clc; clear; close all;  % Clear command window, workspace, and close all figures
addpath("./utils/");      % Add the utils directory to the search path
addpath("./group_analysis/");      % Add the directory to the search path
addpath("./fnc_computation/");      % Add the directory to the search path
addpath("./null_hypothesis/");      % Add the directory to the search path
%% Data Initialization
% Define the number of subjects, timepoints, components, and sampling time
num_subjects = 100;       % Number of subjects
num_timepoints = 150;   % Number of timepoints
num_components = 10;    % Number of components
num_features = nchoosek(num_components, 2); %Unique number of connection pairs
TR = 2;                 % repetition time
bandwidth = [0.01 0.1];
% Generate example post-processed fMRI timecourses
subTcs = randn(num_subjects, num_timepoints, num_components);
subTcs = post_processing(subTcs, TR, bandwidth);

%% Compute Correlation FNC
corr_fnc = compute_correlation_fnc(subTcs);

%% Compute DTW across subjects
win_size = calculate_window_size(TR, bandwidth);
[dtw_stand, dtw_norm] = compute_DTW(subTcs, win_size);

%% Convert DTW distances to FNC
dtw_stand_fnc = -zscore(dtw_stand, 0, "all");
dtw_norm_fnc = -zscore(dtw_norm, 0, "all");

%% Surrogate parameters
num_surrogate = 10;
%% TS surrogate creation
min_shift = round(0.2*num_timepoints); %Set minimum shift to 20% of signal length
subTcs_ts = ts_null_model(subTcs, num_surrogate, min_shift);

%% PR surrogate creation
subTcs_pr = pr_null_model(subTcs, num_surrogate);

%% Computing FNC on surrogate data
surr_corr_fnc = zeros(num_surrogate, num_subjects, num_features);
surr_dtw_stand_fnc = zeros(num_surrogate, num_subjects, num_features);
surr_dtw_norm_fnc = zeros(num_surrogate, num_subjects, num_features);

for surr_num = 1 : num_surrogate
    %Extract single surrogate sample
    surr_subTcs = squeeze(subTcs_pr(:, :, :, surr_num));
    %Perform correlation FNC
    surr_corr_fnc(surr_num, :, :) = compute_correlation_fnc(surr_subTcs);
    %Perform dynamic time warping
    [surr_dtw_stand_dist, surr_dtw_norm_dist] = compute_DTW(surr_subTcs, win_size);
    %For more rigorous comparison, convert each null to fnc's. Can also
    %convert all dtw null into FNC all together for a less rigorous comparison.
    surr_dtw_stand_fnc(surr_num, :, :) = -zscore(surr_dtw_stand_dist, 0, "all");
    surr_dtw_norm_fnc(surr_num, :, :) = -zscore(surr_dtw_norm_dist, 0, "all");
end
%% Perform test between original and surrogates
%Concatenate null standard DTW
surr_dtw_stand_all = reshape(surr_dtw_stand_fnc, num_subjects*num_surrogate, num_features);
%Concatenate null normalized DTW
surr_dtw_norm_all = reshape(surr_dtw_norm_fnc, num_subjects*num_surrogate, num_features);
%Concatenate null correlation
surr_corr_all = reshape(surr_corr_fnc, num_subjects*num_surrogate, num_features);

% null_stats   - A 4 x num_features matrix where:
%              - Row 1 contains uncorrected p-values.
%              - Row 2 contains corrected p-values using Benjamini-
%                Hochberg FDR correction.
%              - Row 3 contains test statistics (z-values from 
%                rank-sum test).
%              - Row 4 contains Cohen's d effect sizes.

dtw_stand_null_stats = perform_null_statistical_test(surr_dtw_stand_all, dtw_stand_fnc);
dtw_norm_null_stats = perform_null_statistical_test(surr_dtw_norm_all, dtw_norm_fnc);
corr_null_stats = perform_null_statistical_test(surr_corr_all, corr_fnc);

%%