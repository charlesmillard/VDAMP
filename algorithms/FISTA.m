function [x_hat, hist] = FISTA(dcoil, mask, x0, opts)
%    The FISTA, S-FISTA (SISTA) or SUREIT algorithm with for
%    image reconstruction from Fourier coefficients, where
%    a sparse model on Haar wavelet coefficients is assumed
%    IN:
%         dcoil: (nx*ny) measurement data i.e. noisy undersampled Fourier coefficients...
%             unsampled coefficients have entry zero
%         mask: (nx*ny) matrix of zeros and ones that signifies the sampling location
%         x0: (nx*ny) ground truth/fully sampled image for reconsturction
%         quality metrics
%         opts: options object with attributes
%             maxIter: maximum number of allowed iterations; default 200
%             verbose: if 1, print progress metrics; default 0
%             scales: number of wavelet decomposition scales; default 4    
%             fistaFlag: 0 for ISTA, 1 for (n-1)/(n+2) version, 2 for
%             (1+sqrt(1+4tk^2))/2 version; default 0
%             saveHist: if 1, save detailed history of reconstruction; default 0
%             SURE: if 1, use SURE to automatically tune thresholds;
%             default 1
%             lambda: if not using SURE, specify sparse weighting
%             sista_div: sista's subband-wise weighting, default 1 (
%
%     OUT:
%         x_hat: VDAMP's x0 reconstruction
%         hist:  history object formed with attributes:
%               timer: (maxIter) time at every iteration             
%               x_mse: (maxIter) ground truth mean-squared error of x              
%               **the following are saved only if saveHist == 1**:
%                   r1: (nx*ny*maxIter) image representation of vector
%                   subject to thresholding (if saveHist == 0, saves last
%                   only)
%                   x: (nx*ny*maxIter) image estimate
%
% The Software does not have 510K approval,  U.S. Food and Drug Administration (FDA) clearance,
% European Medicines Agency (EMA) clearance or CE mark. The Software is intended research purposes only,
% not for diagnostic or clinical use.
%
% Copyright (C) 2019  Charles Millard
% Copyright (C) 2019  Siemens Healthineers
    
    [nx, ny] = size(dcoil);
    
    if (isfield(opts,'maxIter') && ~isempty(opts.maxIter))
        maxIter = opts.maxIter;
    else
        maxIter = 200;
    end
    
    if (isfield(opts,'maxTime') && ~isempty(opts.maxTime))
        maxTime = opts.maxTime;
    else
        maxTime = 60;
    end
    
    if (isfield(opts,'verbose') && ~isempty(opts.verbose))
        verbose = opts.verbose;
    else
        verbose = 0;
    end
    
    if (isfield(opts,'scales') && ~isempty(opts.scales))
        scales = opts.scales;
    else
        scales = 4;
    end
    
    if (isfield(opts,'fistaFlag') && ~isempty(opts.fistaFlag))
        fistaFlag = opts.fistaFlag;
    else
        fistaFlag = 0;
    end
    
    if (isfield(opts,'saveHist') && ~isempty(opts.saveHist))
        saveHist = opts.saveHist;
    else
        saveHist = 1;
    end
    
    if (isfield(opts,'SURE') && ~isempty(opts.SURE)) 
        sureFlag = opts.SURE;
    else
        sureFlag = 1;
    end
    
    if (isfield(opts,'sista_div') && ~isempty(opts.sista_div) && iscell(opts.sista_div)) 
        sista_div = opts.sista_div;
    else
        sista_div = 1;
    end
    
    if (isfield(opts,'wavType') && ~isempty(opts.wavType)) 
        wavType = opts.wavType;
    else
        wavType = 'haar';
    end
    
    
    if sureFlag  == 0 
         if (isfield(opts,'lambda') && ~isempty(opts.lambda))
             sparse_weight = opts.lambda;
         else
             error('SURE is off and no lambda has been selected: please set opts.lambda')
         end
    end
    

    r2 = cell(scales, 1);
    lambda = cell(scales, 1);
    tau = cell(scales, 1); 
    for s =1:scales
        r2{s} = struct('H', [], 'V', [], 'D', []);   
        if sureFlag == 0
            if iscell(sista_div)
                lambda{s} = struct('H', sparse_weight/(2*sista_div{s}.H), 'V', sparse_weight/(2*sista_div{s}.V), 'D', sparse_weight/(2*sista_div{s}.D));
            else
                lambda{s} = struct('H', sparse_weight, 'V', sparse_weight, 'D', sparse_weight);
            end
        end
    end
    
    r2{scales}.A = [];
    if sureFlag == 0
        if iscell(sista_div)
            lambda{scales}.A = sparse_weight./(2*sista_div{scales}.A);
        else
            lambda{scales}.A = sparse_weight;
        end
    end
    
    W0 = multiscaleDecomp(x0, scales, wavType);

    r = ifftnc(mask.*dcoil);
    
    if iscell(sista_div)
        z_wav = multiscaleDecomp(r, scales, wavType);
        for s=1:scales
            subbands = fieldnames(z_wav{s});           
                for i = 1:numel(subbands)
                    band_name = subbands{i};
                    C{s}.(band_name) = z_wav{s}.(band_name)./sista_div{s}.(band_name);
                end           
        end
    else    
        C = multiscaleDecomp(r, scales, wavType);
    end
    
    var_est = immse(pyramid(C), pyramid(W0));
    
    for s =1:scales
        tau{s} = struct('H', var_est, 'V', var_est, 'D', var_est);
    end
    tau{scales}.A = var_est;
    
    if saveHist
        hist.r1 = zeros(nx,ny,maxIter);
        hist.x = zeros(nx,ny,maxIter);
        hist.timer = zeros(1, maxIter);
    end
    
    t_new = 1; % fista weighting initialisation
    x = 0;
    hist.x_mse = zeros(1, maxIter);
    
    time_init = tic; 
    
    for iter = 1:maxIter
        
        if saveHist
            hist.r1(:,:,iter) = r;
        end
        
        if verbose
            disp(['                 ITERATION ', num2str(iter)])
            fprintf( ...
                'True RMSE = %f\n', ...
                sqrt(mean(abs(r(:)-x0(:)).^2)));
        end
        
        % thresholding
        if sureFlag
            C_thr = multiscaleSUREsoft(C, tau);
        else      
            C_thr = multiscaleComplexSoft(C, tau, lambda);
        end
        
        x_new = multiscaleRecon(C_thr, wavType);
        x_old = x;
      
        if fistaFlag == 1
            a = 1+ (iter-1)/(iter+2);
            b = (1 - iter)/(iter+2);
            x = a*x_new + b*x_old;
        elseif fistaFlag == 2
            t_old = t_new;
            t_new = (1+sqrt(1+4*t_old^2))/2;
            x = x_new + (t_old -1)/t_new*(x_new -x_old);
        else    
            x = x_new;
        end
        
        % gradient step
        z = dcoil - mask.*fftnc(x);       
        if iscell(sista_div)
            z_wav = multiscaleDecomp(ifftnc(z), scales, wavType);
            x_wav = multiscaleDecomp(x, scales, wavType);
            for s=1:scales
                subbands = fieldnames(C_thr{s});           
                    for i = 1:numel(subbands)
                        band_name = subbands{i};
                        
                        C{s}.(band_name) = z_wav{s}.(band_name)./sista_div{s}.(band_name) + x_wav{s}.(band_name);
                    end
            end
        else    
            r = x+ifftnc(z);
            C = multiscaleDecomp(r, scales, wavType);
        end
        
        x_hat = ifftnc(mask.*z) + x_new; % current estimate
        hist.x_mse(iter) = immse(x0, x_hat);
 
        hist.timer(iter) = toc(time_init); 
        
        if saveHist
            hist.x(:,:,iter) = r;
        end
    
        if hist.timer(iter)> maxTime
            break
        end
        
        var_est =  immse(pyramid(C), pyramid(W0)); 
        for s =1:scales
            tau{s} = struct('H', var_est, 'V', var_est, 'D', var_est);
        end
        tau{scales}.A = var_est;  
    end

    if ~saveHist
        hist.r1 = r;
    end
end

