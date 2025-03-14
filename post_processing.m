function [pp_subTcs, pp_subTcs_nz] = post_processing(subTcs, TR, band)
    % post_processing performs quality control measures on the fMRI data;
    % detrending; despiking; filtering; zscoring
    %  
    % Inputs:
    % 1. subTcs    - The fMRI time series across subjects/individuals (Subject by Time by Component)
    % 2. TR        - sampling rate of the time series (i.e., Tr in fMRI)
    % 3. band      - frequency band for filtering (i.e., 0.01hz to 0.15hz for BOLD signal)
    %
    % Outputs:
    % 1. pp_subTcs - detrended, despiked, filtered, and zscored fMRI time
    % series (Subject by Time by Component)
    % 2. pp_subTcs_nz - detrended, despiked, and filtered, NO z-score!! (Subject by Time by Component)
    %  
    % Examples: 
    %   Here 10 time series are generated by drawing samples from a normal
    %   distribution for 10 Subjects. Asumming a TR of 1 second and filter band of 0.01hz to 0.15hz 
    %  
    %   subTcs = randn(10,1000,10);
    %   TR = 1;
    %   band = [0.01 0.15];
    %   [pp_subTcs, pp_subTcs_nz] = post_processing(subTcs,TR,band);

    pp_subTcs = zeros(size(subTcs, 1), size(subTcs, 2), size(subTcs, 3));
    pp_subTcs_nz = zeros(size(subTcs, 1), size(subTcs, 2), size(subTcs, 3));
    srate = 1/TR;
    
    %For each subject, perform all four steps of quality control
    % measures
    for sub = 1:size(subTcs, 1)
        fprintf('Processing subject: %d\n', sub)
        pp_subTcs_nz(sub, :, :) = detrend(squeeze(subTcs(sub, :, :))); % Linear Detrending

        pp_subTcs_nz(sub, :, :) = filtering(squeeze(pp_subTcs_nz(sub, :, :)), srate, band); % Filtering
        
        pp_subTcs(sub, :, :) = zscore(squeeze(pp_subTcs_nz(sub, :, :))); % Z-scoring
    end  
end

function filtered_Tc = filtering(subTc, srate, band)
    nyquist = srate/2;
    Wp = band/nyquist;
    Ws = [band(1)*0.5 band(2)*1.5]/nyquist;
    Rp = 3;
    Rs = 30;
    [n,Wn] = buttord(Wp,Ws,Rp,Rs);

    % filter coefficients
    [fkernB,fkernA] = butter(n, Wn);

    % apply the filter through all components
    subTc_zeros = zeros(size(subTc, 1), size(subTc, 2));
    for sub = 1 : size(subTc, 2)
        subTc_zeros(:,sub) = filtfilt(fkernB,fkernA,subTc(:,sub));
    end
    filtered_Tc = subTc_zeros;
end