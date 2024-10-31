function d = calculate_cohen_d(x, y)
% CALCULATE_COHEN_D Compute Cohen's d effect size between two samples.
%
%   d = CALCULATE_COHEN_D(x, y) calculates Cohen's d, a measure of effect size, 
%   for two independent samples. Cohen's d is the difference in means 
%   between two samples divided by the pooled standard deviation, 
%   giving a standardized effect size.
%
%   INPUTS:
%       x - A vector of sample data for group 1.
%       y - A vector of sample data for group 2.
%
%   OUTPUT:
%       d - Cohen's d effect size, representing the standardized difference 
%           between the two groups.
%
%   Example:
%       % Calculate Cohen's d for two sample groups
%       group1 = randn(50, 1);
%       group2 = randn(50, 1) + 0.5; % Mean shifted by 0.5 for illustration
%       d = calculate_cohen_d(group1, group2);
%
%   Notes:
%       - Cohen's d is commonly used to assess effect size in hypothesis testing.
%       - Assumes equal variances between the two groups (pooled standard deviation).
%
%   Author: [Sir-Lord]
%   Version: 1.0

    % Sample sizes for each group
    n1 = length(x);
    n2 = length(y);

    % Means for each group
    mean_x = mean(x);
    mean_y = mean(y);

    % Standard deviations for each group
    std_x = std(x);
    std_y = std(y);

    % Pooled standard deviation calculation
    pooled_std = sqrt(((n1 - 1) * std_x^2 + (n2 - 1) * std_y^2) / (n1 + n2 - 2));
    
    % Cohen's d calculation
    d = (mean_x - mean_y) / pooled_std;
end
