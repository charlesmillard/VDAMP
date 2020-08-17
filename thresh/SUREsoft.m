function [gz, df] = SUREsoft(z, V)
    % Complex soft thresholding with a threshold chosen by SURE
    % IN:
    %   z: vector corrupted with white complex Gaussian noise
    %   V: variance of z
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

    
    z2 = z(:);
    [mag, index] = sort(abs(z2), 'descend');
    
    V = V(:) .*ones(size(z2));
    V = V(index);
    
    z0 = ones(size(mag));
    
    lambda = mag;
    
    SURE_inf = flipud(cumsum(flipud(mag.^2)));
    SURE_sup = cumsum(z0).*lambda.^2 - lambda.*cumsum(V./mag) + 2*cumsum(V);
    SURE = SURE_inf + SURE_sup - sum(V);
    
    [val, idx] = min(SURE);
    [gz, df] = complexSoft(z, lambda(idx));
end

