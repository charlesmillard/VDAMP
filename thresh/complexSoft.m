function [gz, df] = complexSoft(z, lambda)
    % Complex soft thresholding operator with threshold lambda
    % IN:
    %   z: noisy vector
    %   lambda: threshold
    % OUT:
    %   gz: thresholded z
    %   df: number of degrees of freedom
    %
    % The Software does not have 510K approval,  U.S. Food and Drug Administration (FDA) clearance,
    % European Medicines Agency (EMA) clearance or CE mark. The Software is intended research purposes only,
    % not for diagnostic or clinical use.
    %
    % Copyright (C) 2019  Charles Millard
    % Copyright (C) 2019  Siemens Healthineers
    
    mag = abs(z);
    gdual = min(lambda ./ mag, 1);
    gz = z.*(1 - gdual);
    
    % To compute the degrees of freedom of a complex function, consider it 
    % as a function from R^2 to R^2:
    %
    % df = d(real(gz))/d(real(z)) + d(imag(gz))/d(imag(z))
    %
    
    df = 2 - (2 - (gdual<1)).*gdual;
end

