function [q_values, p_values] = compute_glm_group_difference(fnc_matrix, group, covariates)
% compute_glm_group_difference - Computes group differences across FNC features using GLM.
%
% Syntax:
%   [q_values, p_values] = compute_glm_group_difference(fnc_matrix, group, covariates)
%
% Inputs:
%   fnc_matrix - A matrix (n x m), where n is the number of observations (e.g., subjects)
%                and m is the number of FNC features. Each column represents an FNC feature.
%
%   group      - A binary vector (n x 1) representing the group assignment of each subject.
%                Typically, values are 0 and 1, where 1 might represent a specific group (e.g., patients)
%                and 0 another (e.g., controls).
%
%   covariates - A table (n x p) of covariates for each subject, where n is the number of
%                observations and p is the number of covariates. This table can include any
%                relevant variables, such as age, sex, site, or other covariates.
%
% Outputs:
%   q_values   - A vector (m x 1) of q-values for each FNC feature after applying
%                Benjamini-Hochberg FDR correction.
%
%   p_values   - A vector (m x 1) of p-values for each FNC feature before FDR correction.
%
% Example:
%   % Example Data (Randomized for demonstration)
%   n = 311;          % Number of subjects
%   m = 100;          % Number of FNC features
%   fnc_matrix = randn(n, m);       % Example FNC matrix (311 subjects, 100 features)
%   group = randi([0, 1], n, 1);    % Random binary group assignment
%   
%   % Creating a table of covariates
%   age = randi([20, 80], n, 1);    % Random ages between 20 and 80
%   sex = randi([0, 1], n, 1);      % Random binary sex assignment
%   site = randi([1, 3], n, 1);     % Random site assignment with 3 possible sites
%   covariates = table(age, sex, site); % Table of covariates
%   
%   % Call the function
%   [q_values, p_values] = compute_glm_group_difference(fnc_matrix, group, covariates);

    % Validate inputs
    if height(covariates) ~= size(fnc_matrix, 1)
        error('The number of rows in covariates must match the number of rows in fnc_matrix.');
    end
    
    % Determine the number of FNC features
    num_features = size(fnc_matrix, 2);
    
    % Initialize p-values vector for each feature
    p_values = zeros(num_features, 1);
    
    % Add 'group' as the first column in the covariates table
    covariates = addvars(covariates, group, 'NewVariableNames', 'group', 'Before', 1);
    
    % Loop through each FNC feature and fit a GLM
    for feature_idx = 1:num_features
        % Extract the current FNC feature data (a column vector)
        fnc_feature = fnc_matrix(:, feature_idx);
        
        % Add FNC feature as the response variable to the covariates table
        data = addvars(covariates, fnc_feature, 'NewVariableNames', 'fnc');
        
        % Fit a GLM with 'group' as the predictor and covariates as control variables
        mdl = fitglm(data, 'ResponseVar', 'fnc');
        
        % Extract the p-value for the 'group' effect (assumed to be the 2nd coefficient)
        p_values(feature_idx) = mdl.Coefficients.pValue('group');
    end
    
    % Apply Benjamini-Hochberg FDR correction
    q_values = mafdr(p_values, 'BHFDR', true);
end
