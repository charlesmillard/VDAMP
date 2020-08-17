function bands = multiscaleDecomp(A, scales)
% Perform a Haar wavelet decomposition    
% IN: 
%   A: image
%   scales: number of decomposition scales
% OUT:
%   bands: cell array of wavelet bands
%
% The Software does not have 510K approval,  U.S. Food and Drug Administration (FDA) clearance,
% European Medicines Agency (EMA) clearance or CE mark. The Software is intended research purposes only,
% not for diagnostic or clinical use.
%
% Copyright (C) 2019  Charles Millard
% Copyright (C) 2019  Siemens Healthineers

    bands = cell(scales, 1);
    
    for s=1:scales
        [A, H, V, D] = dwt2(A, 'haar');
        bands{s} = struct( ...
            'H', H, ...
            'V', V, ...
            'D', D);
    end
    
    bands{scales}.A = A;
end