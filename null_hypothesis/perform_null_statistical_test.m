function null_stats = perform_null_statistical_test(orig_data, surr_data)
% PERFORM_NULL_STATISTICAL_TEST Conduct statistical testing on null data.
%
%   null_stats = PERFORM_NULL_STATISTICAL_TEST(orig_data, surr_data)
%   computes statistical metrics for each feature by comparing surrogate
%   data against original data. The function outputs corrected and uncorrected 
%   p-values, test statistics, and Cohen's d effect sizes for 
%   each feature.
%
%   INPUTS:
%       orig_data    - A matrix or array of original observed data with 
%                      dimensions [observations x features].      
%       surr_data    - A matrix or array of surrogate data with dimensions 
%                      [observations x features].
%
%       
%
%   OUTPUT:
%       null_stats   - A 4 x num_features matrix where:
%                      - Row 1 contains uncorrected p-values.
%                      - Row 2 contains corrected p-values using Benjamini-
%                        Hochberg FDR correction.
%                      - Row 3 contains test statistics (z-values from 
%                        rank-sum test).
%                      - Row 4 contains Cohen's d effect sizes.
%
%   Example:
%       % Conduct a null statistical test with 100 features
%       null_stats = perform_null_statistical_test(orig_data, surr_data);
%
%   Notes:
%       - Performs rank-sum test for each feature comparison between 
%         surrogate and original data.
%       - Corrects p-values using the Benjamini-Hochberg FDR procedure.
%
%   Author: [Sir-Lord]
%   Version: 1.0

    % Initialize output matrices and variables
    num_features = size(orig_data, 2);
    null_stats = zeros(4, num_features);          % To store all statistics
    uncorrected_p_values = zeros(1, num_features);
    t_values = zeros(1, num_features);
    cohens_d_effect_size = zeros(1, num_features);

    % Loop through each feature to perform statistical testing
    for f = 1:num_features
        % Extract data for the current feature
        x = orig_data(:, f);
        y = surr_data(:, f); 

        % Conduct rank-sum test and extract p-value and z-value (test statistic)
        [p_val, ~, stats] = ranksum(x,y);
        t_values(f) = stats.zval;

        % Store uncorrected p-value and calculate Cohen's d effect size
        uncorrected_p_values(f) = p_val;
        cohens_d_effect_size(f) = calculate_cohen_d(x, y);
    end
    
    % Apply Benjamini-Hochberg FDR correction to p-values
    corrected_p_values = mafdr(uncorrected_p_values, "BHFDR", true);
    
    % Populate output matrix with calculated statistics
    null_stats(1, :) = corrected_p_values;
    null_stats(2, :) = uncorrected_p_values;
    null_stats(3, :) = t_values;
    null_stats(4, :) = cohens_d_effect_size;

end