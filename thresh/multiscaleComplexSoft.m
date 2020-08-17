function [bands, err, df] = multiscaleComplexSoft(bands, var, lambda)
    % Soft thresholding of a wavelet respresentation with fixed sparse
    % weighting
    % IN:
    %   bands: wavelet coefficients
    %   var: estimated MSE of bands
    %   lambda: sparse weighting
    % OUT:
    %   bands: thresholded wavelet coefficients
    %   err: estimated MSE of bands
    %   number of degrees of freedom of bands
    %
    % The Software does not have 510K approval,  U.S. Food and Drug Administration (FDA) clearance,
    % European Medicines Agency (EMA) clearance or CE mark. The Software is intended research purposes only,
    % not for diagnostic or clinical use.
    %
    % Copyright (C) 2019  Charles Millard
    % Copyright (C) 2019  Siemens Healthineers
    
    scales = numel(bands);
    df = bands;
    err = cell(scales, 1);
    
    for s = 1:scales
        subbands = fieldnames(bands{s});
        
        for i = 1:numel(subbands)
            band_name = subbands{i};
            
            [bands{s}.(band_name), df{s}.(band_name)] = ...
                complexSoft(bands{s}.(band_name), ...
                var{s}.(band_name) .* lambda{s}.(band_name));
            
            e = df{s}.(band_name) .* var{s}.(band_name);
            
            % Dive by 2 because the input variance is the expected energy 
            % of the complex noise, i.e. sqrt(2) higher than the variance
            % of the real and imaginary parts
            err{s}.(band_name) = mean(e(:)) / 2;
        end
    end
end

