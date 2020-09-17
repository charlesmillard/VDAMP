function [meank_re, meank_im] = subband_kurtosis(w0, r1, wavType)
% calculates the real and imaginary parts of the mean subband-wise kurtosis of r1
%
%
% The Software does not have 510K approval,  U.S. Food and Drug Administration (FDA) clearance,
% European Medicines Agency (EMA) clearance or CE mark. The Software is intended for research purposes only,
% not for diagnostic or clinical use.
%
% Copyright (C) 2020  Charles Millard
% Copyright (C) 2020  Siemens Healthineers

scales = numel(w0);

k = zeros(3,scales);
k_im = zeros(3,scales);



C = multiscaleDecomp(r1,scales, wavType);

for s = 1:scales
    k(1,s) =  kurtosis(real(C{s}.H(:) - w0{s}.H(:)));
    k(2,s) =  kurtosis(real(C{s}.V(:) - w0{s}.V(:)));
    k(3,s) =  kurtosis(real(C{s}.D(:) - w0{s}.D(:)));

    k_im(1,s) =  kurtosis(imag(C{s}.H(:) - w0{s}.H(:)));
    k_im(2,s) =  kurtosis(imag(C{s}.V(:) - w0{s}.V(:)));
    k_im(3,s) =  kurtosis(imag(C{s}.D(:) - w0{s}.D(:)));
end
kA =  kurtosis(real(C{s}.A(:) - w0{s}.A(:)));
kA_im =  kurtosis(imag(C{s}.A(:) - w0{s}.A(:)));

meank_re = mean([k(:); kA] -3); % subtract 3 for excess kurtosis
meank_im = mean([k_im(:); kA_im] - 3);

% disp(['Mean kurtosis: real: ', num2str(meank_re), ' || imag: ', num2str(mean(meank_im))])