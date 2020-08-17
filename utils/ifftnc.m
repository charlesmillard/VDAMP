function x = ifftnc(k)
% inverse FFT of a k with low freq Fourier coefficients at centre, normalised 
% IN:
%   k: Fourier coefficients
% OUT:
%   x: inverse FFT of k
%
% The Software does not have 510K approval,  U.S. Food and Drug Administration (FDA) clearance,
% European Medicines Agency (EMA) clearance or CE mark. The Software is intended research purposes only,
% not for diagnostic or clinical use.
%
% Copyright (C) 2019  Charles Millard
% Copyright (C) 2019  Siemens Healthineers

x = ifftshift(ifftn(ifftshift(k)));
x = sqrt(length(k(:)))*x;
