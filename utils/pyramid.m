function img = pyramid(bands)
% full wavelet representation 
% IN:
%   bands: wavelet coefficients as cell array 
% OUT:
%   img: pyramidal representation of wavelet   
%
% The Software does not have 510K approval,  U.S. Food and Drug Administration (FDA) clearance,
% European Medicines Agency (EMA) clearance or CE mark. The Software is intended research purposes only,
% not for diagnostic or clinical use.
%
% Copyright (C) 2019  Charles Millard
% Copyright (C) 2019  Siemens Healthineers

    scales = numel(bands);
    img = bands{scales}.A;
    
    for s = scales:-1:1
        img = [img bands{s}.H; bands{s}.V bands{s}.D];
    end
end

