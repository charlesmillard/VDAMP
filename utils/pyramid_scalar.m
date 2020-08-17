function img = pyramid_scalar(bands, nx, ny)
% converts a cell with one scalar per subband to the full wavlet representation
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
    img = bands{scales}.A.*ones(nx*2^-(scales),ny*2^-(scales));
    
    for s = scales:-1:1
        bandH = bands{s}.H.*ones(nx*2^-(s),ny*2^-(s));
        bandV = bands{s}.V.*ones(nx*2^-(s),ny*2^-(s));
        bandD = bands{s}.D.*ones(nx*2^-(s),ny*2^-(s));
        img = [img, bandH ; bandV,  bandD];
    end
end

