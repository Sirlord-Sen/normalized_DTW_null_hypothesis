function [dtw_stand_dist, dtw_norm_dist] = compute_DTW(subTcs, win_size)
    % compute_DTW - Calculates standard and normalized DTW distances between
    % components for multiple subjects and stores the results in vectors.
    %
    % Syntax:
    %   [dtw_stand_dist, dtw_norm_dist] = compute_DTW(subTcs, win_size)
    %
    % Inputs:
    %   subTcs   - 3D matrix of subject timecourses with dimensions:
    %              [num_subjects, num_timepoints, num_components].
    %   win_size - Window size for DTW. Must be a positive integer.
    %
    % Outputs:
    %   dtw_stand_dist - Vectorized standard DTW distances across components.
    %   dtw_norm_dist  - Vectorized normalized DTW distances across components.
    %
    % Example:
    %   [dtw_stand_dist, dtw_norm_dist] = compute_DTW(subTcs, 5);
    
    % Input validation
    if nargin ~= 2
        error('compute_DTW requires exactly two input arguments: subTcs and win_size.');
    end
    
    % Check that subTcs is a 3D numeric array
    if ~isnumeric(subTcs) || ndims(subTcs) ~= 3
        error('subTcs must be a 3D numeric array with dimensions [num_subjects, num_timepoints, num_components].');
    end

    % Check that win_size is a positive integer
    if ~isscalar(win_size) || ~isnumeric(win_size) || win_size <= 0 || mod(win_size, 1) ~= 0
        error('win_size must be a positive integer.');
    end

    % Extract dimensions
    [num_subjects, num_timepoints, num_components] = size(subTcs);

    % Initialize matrices for DTW distances
    dtw_stand = zeros(num_subjects, num_components, num_components);
    dtw_norm = zeros(num_subjects, num_components, num_components);

    % Compute DTW distances
    for sub_num = 1:num_subjects
        fprintf('DTW subject: %d of %d\n', sub_num, num_subjects);

        % Only fill the lower triangular part
        for comp1 = 2:num_components
            for comp2 = 1:comp1-1
                % Extract and z-score the timecourses for each component pair
                x = zscore(squeeze(subTcs(sub_num, :, comp1)));
                y = zscore(squeeze(subTcs(sub_num, :, comp2)));

                % Calculate DTW distance and obtain the matching indices
                [dtw_d, ix, iy] = dtw(x, y, win_size);

                % Store standard and normalized DTW distances
                dtw_stand(sub_num, comp1, comp2) = dtw_d;
                dtw_norm(sub_num, comp1, comp2) = dtw_d / length(ix);
            end
        end
    end

    % Vectorize the DTW matrices using the icatb_mat2vec function
    dtw_stand_dist = icatb_mat2vec(dtw_stand);
    dtw_norm_dist = icatb_mat2vec(dtw_norm);
end
