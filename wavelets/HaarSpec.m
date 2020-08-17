function spec = HaarSpec(scales, sample)
    % Compute power spectrum of Haar wavelet filters at multiple scales
    % IN:
    %   scales: number of decomposition scales
    %   sample: sample dimensionality
    % OUT:
    %   spec: power spectrum
    %
    % The Software does not have 510K approval,  U.S. Food and Drug Administration (FDA) clearance,
    % European Medicines Agency (EMA) clearance or CE mark. The Software is intended research purposes only,
    % not for diagnostic or clinical use.
    %
    % Copyright (C) 2019  Charles Millard
    % Copyright (C) 2019  Siemens Healthineers
    
    spec = zeros(sample, scales, 2);
    
    filtLen = 1;
    
    for s=1:scales
        % compute band-pass
        filtLen = filtLen*2;
        
        filt = zeros(sample, 1);
        filt(1:filtLen) = 2^(-s/2);
        
        spec(:, s, 1) = abs(fftnc(filt)).^2;
        
        filt = zeros(sample, 1);
        filt(1:filtLen/2) = -2^(-s/2);
        filt(filtLen/2+1:filtLen) = 2^(-s/2);
        
        spec(:, s, 2) = abs(fftnc(filt)).^2;
    end
end