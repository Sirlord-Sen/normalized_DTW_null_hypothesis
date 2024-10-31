function surr_data = ts_null_model(data, num_surrogate, min_shift)
% TS_NULL_MODEL Generate temporally shifted surrogate time series.
%
%   surr_data = TS_NULL_MODEL(data, num_surrogate, max_shift) generates 
%   surrogate data by performing temporal shifts on the input time courses. 
%   The number of surrogates to generate and the maximum shift allowed 
%   are specified as inputs.
%
%   INPUTS:
%       data         - A 3D numeric array of input time courses with 
%                      dimensions [subjects x time points x components].
%
%       num_surrogate - The number of surrogate data sets to generate.
%
%       min_shift    - The minimun number of time points by which each 
%                      component is allowed to shift. The shift is randomly
%                      selected within the range [min_shift, time points-min_shift].
%
%   OUTPUTS:
%       surr_data    - A 4D numeric array of surrogate data with 
%                      dimensions [subjects x time points x components x num_surrogate].
%
%   Example:
%       % Generate surrogate data for 10 subjects, 1000 time points, 
%       % 3 components, 50 surrogates, and a minimum shift of 20% of time points.
%       data = randn(10, 1000, 3);   % Example input data
%       num_surrogate = 50;          
%       min_shift = round(0.2*size(data, 2));             
%       surr_data = ts_null_model(data, num_surrogate, max_shift);
%
%   Notes:
%       - Temporal shifting is performed independently for each subject 
%         and component.
%       - Shift amounts are randomly generated within the specified range 
%          [min_shift, time points-min_shift].
%
%   Author: [Sir-Lord]
%   Version: 1.0

    % Input validation
    if ~isnumeric(data) || ndims(data) ~= 3
        error('Input "data" must be a 3-dimensional numeric array.');
    end

    if ~isnumeric(num_surrogate) || ~isscalar(num_surrogate) || num_surrogate <= 0
        error('Input "num_surrogate" must be a positive scalar integer.');
    end

    if ~isnumeric(min_shift) || ~isscalar(min_shift) || min_shift <= 0
        error('Input "min_shift" must be a positive scalar integer.');
    end

    % Extract dimensions of the input matrix
    [subs, tps, comps] = size(data);

    % Initialize the output matrix with zeros
    surr_data = zeros(subs, tps, comps, num_surrogate);
    
    % Generate the surrogates with temporal shifting
    for surr_num = 1:num_surrogate
        fprintf("Generating surrogate: %d/%d\n", surr_num, num_surrogate);
        
        % Loop through each subject
        for sub = 1:subs
            % Loop through each component
            for comp = 1:comps
                % Randomly determine shift amount within range
                shift_amount = randi([min_shift, tps-min_shift]);
                
                % Perform circular shift on the data
                surr_data(sub, :, comp, surr_num) = circshift(squeeze(data(sub, :, comp)), shift_amount);    
            end 
        end
    end
end
