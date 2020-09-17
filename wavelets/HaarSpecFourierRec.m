function spec = WavSpec(scales, sample, filt)
    % Compute power spectrum of Haar wavelet filters at multiple scales
    
    if nargin < 3
        filt = [1, -1; 1, 1] / sqrt(2);
    end
    
    spec = zeros(sample, scales, 2);
    
    L = zeros(sample, 1);
    H = zeros(sample, 1);
    
    L(1:size(filt, 1)) = filt(:, 1);
    H(1:size(filt, 1)) = filt(:, 2);
    
    spec(:,1,1) = abs(fft(L)).^2;
    spec(:,1,2) = abs(fft(H)).^2;
    
    sigLen = sample;
    numBlock = 1;
    
    for s=2:scales
        sigLen = sigLen / 2;
        numBlock = numBlock *2;
        
        % Multiply with the previous low-pass spectrum with aliasing
        spec(:, s, 1) = spec(:,s-1,1) .* ...
            reshape(spec(1:numBlock:end,1,1)*ones(1,numBlock), [sample 1]);
        spec(:, s, 2) = spec(:,s-1,1) .* ...
            reshape(spec(1:numBlock:end,1,2)*ones(1,numBlock), [sample 1]);
    end
    
    spec = fftshift(spec, 1) / sample;
end