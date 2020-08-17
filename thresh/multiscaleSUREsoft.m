function [bands, err, df] = multiscaleSUREsoft(bands, var)
    % Soft thresholding of an entire wavelet respresentation with soft
    % threshold tuned by SURE
    % IN:
    %   bands: wavelet coefficients
    %   var: estimated MSE of bands
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
    err = bands;
    df = bands;
    
    for s = 1:scales
        subbands = fieldnames(bands{s});
        
        for i = 1:numel(subbands)
            band_name = subbands{i};
            
            [bands{s}.(band_name), df{s}.(band_name)] = ...
                SUREsoft(bands{s}.(band_name), ...
                var{s}.(band_name));
            
            e = df{s}.(band_name) .* var{s}.(band_name);
            err{s}.(band_name) = mean(e(:)) / 2;
        end
    end
end

