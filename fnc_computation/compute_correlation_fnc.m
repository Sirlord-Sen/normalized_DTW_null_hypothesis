function [corr_fnc] = compute_correlation_fnc(subTcs)
    % compute_correlation_fnc - Calculates correlation between
    % components for multiple subjects and stores the results in a vector.
    %
    % Syntax:
    %   [corr_fnc] = compute_correlation_fnc(subTcs)
    %
    % Inputs:
    %   subTcs   - 3D matrix of subject timecourses with dimensions:
    %              [num_subjects, num_timepoints, num_components].
    %
    % Outputs:
    %   corr_fnc - Vectorized standard correlation fnc across components.
    %
    % Example:
    %   [corr_fnc] = compute_correlation_fnc(subTcs);
    
    % Input validation
    if nargin ~= 1
        error('compute_correlation_fnc requires exactly 1 input arguments: subTcs.');
    end
    
    % Check that subTcs is a 3D numeric array
    if ~isnumeric(subTcs) || ndims(subTcs) ~= 3
        error('subTcs must be a 3D numeric array with dimensions [num_subjects, num_timepoints, num_components].');
    end

    % Extract dimensions
    [num_subjects, num_timepoints, num_components] = size(subTcs);
    num_features = nchoosek(num_components, 2);
    % Initialize matrices for DTW distances
    corr_fnc = zeros(num_subjects, num_features);

    % Compute DTW distances
    for sub_num = 1:num_subjects
        fprintf('Correlation subject: %d of %d\n', sub_num, num_subjects);

        corr_fnc(sub_num, :) = icatb_mat2vec(corr(squeeze(subTcs(sub_num, :, :))));
    end

end
