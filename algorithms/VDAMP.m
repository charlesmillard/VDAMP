function [x_hat, hist] = VDAMP(dcoil, mask, prob_map, var0, x0, opts)
%     The Variable Density Approximate Message Passing (VDAMP) algorithm for
%     reconstruction of a natural image from Fourier coefficients, where
%     a sparse model on Haar wavelet coefficients is assumed
%     IN:
%         dcoil: (nx*ny) measurement data i.e. noisy undersampled Fourier coefficients...
%             unsampled coefficients have entry zero
%         mask: (nx*ny) matrix of zeros and ones that signifies of sampling location
%         prob_map: (nx*ny) matrix of sampling probabilities
%         var0: variance of measurement noise 
%         x0: (nx*ny) ground truth/fully sampled image for tracking reconstruction progress
%         opts: options object with attributes
%             maxIter: maximum number of allowed iterations; default 50
%             verbose: if 1, print progress metrics; default 0
%             scales: number of wavelet decomposition scales; default 4    
%             saveHist: if 1, save detailed history of reconstruction; default 0
%             SURE: if 1, use SURE to automatically tune thresholds;
%             default 1
%             lambda: if not using SURE, specify sparse weighting
%             denoiserDiv: if 1, use SURE to select denoising divisor Cdiv,
%             otherwise use (1-alpha) 
%
%      OUT:
%         x_hat: VDAMP's x0 estimate
%         hist:  history object formed with attributes:
%               timer: (maxIter) time at every iteration
%               Cthr: (nx,ny,maxIter) post-thresholding estimate in wavelet domain
%               x_mse: (maxIter) ground truth mean-squared error of x              
%               **the following are saved only if saveHist == 1**:
%                   r1: (nx*ny*maxIter) image representation of vector subject to thresholding
%                   true_err_C: (maxIter*4*scales) band-wise ground truth RMSE of r, order H, V, D, A...
%                            Those scales without A have zeros
%                   belief_std_C: (maxIter*4*scales) estimate of band-wise RMSE of r
%                   RMSEyest: (maxIter) root-mean of tau i.e. estimate of RMSE of r
%                   true_err_Cthr: (maxIter*4*scales) band-wise ground truth RMSE of w (thresholded r1)
%                   belief_std_Cthr: (maxIter*4*scales) estimate of band-wise RMSE of w         
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
        maxIter = 50;
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
    
    if (isfield(opts,'saveHist') && ~isempty(opts.saveHist))
        saveHist = opts.saveHist;
    else
        saveHist = 0;
    end
    
    if (isfield(opts,'SURE') && ~isempty(opts.SURE)) 
        sureFlag = opts.SURE;
    else
        sureFlag = 1;
    end
    
    if (isfield(opts,'denoiserDiv') && ~isempty(opts.denoiserDiv)) 
        denoiserDiv = opts.denoiserDiv;
    else
        denoiserDiv = 1;
    end
    
    if sureFlag  == 0 
         if (isfield(opts,'lambda') && ~isempty(opts.lambda))
             sparse_weight = opts.lambda;
         else
             error('SURE is off and no lambda has been selected: please set opts.lambda')
         end
    end
    
        
    tau = cell(scales, 1); 
    lambda = cell(scales, 1);
    err = cell(scales, 1);
    C_tilde= cell(scales, 1);
    alpha = cell(scales, 1);
    for s =1:scales
        tau{s} = struct('H', 0, 'V', 0, 'D', 0);
        err{s} = struct('H', 0, 'V', 0, 'D', 0);
        C_tilde{s} = struct('H', [], 'V', [], 'D', []);
        alpha{s} = struct('H', 0, 'V', 0, 'D', 0);       
        if sureFlag ==0
            lambda{s} = struct('H', sparse_weight, 'V', sparse_weight, 'D', sparse_weight);
        end       
    end  
    tau{scales}.A = 0;
    err{scales}.A = 0;
    C_tilde{scales}.A = [];
    alpha{s}.A = 0;   
    if sureFlag == 0
        lambda{scales}.A = sparse_weight/100;
    end
    
    specX = HaarSpec(scales, nx);
    specY = HaarSpec(scales, ny);
    
    W0 = multiscaleDecomp(x0, scales);
    
    
    inv_p = prob_map.^(-1);
    inv_p_m1 = inv_p - 1;
    
    % unbaised initialisation
    m2 = inv_p.*mask.*dcoil;
    r = ifftnc(m2);
    C = multiscaleDecomp(r, scales);
    
    % calculate initial tau
    y_cov = mask.*( inv_p_m1 .* inv_p .* abs(dcoil).^2 + inv_p.*var0); % k-space error model
    for s = 1:scales
        tau{s}.H = specX(:, s, 2)'* y_cov * specY(:, s, 1); 
        tau{s}.V = specX(:, s, 1)'* y_cov * specY(:, s, 2); 
        tau{s}.D = specX(:, s, 2)'* y_cov * specY(:, s, 2); 
    end
    tau{scales}.A = specX(:, s, 1)'* y_cov * specY(:, s, 1); 
    
   
    % prepare history
    if saveHist    
        hist.x = zeros(nx,ny,maxIter);
        hist.true_err_C= zeros(maxIter, 4, scales);
        hist.belief_std_C = zeros(maxIter, 4, scales);
        hist.true_err_Cthr= zeros(maxIter, 4, scales);
        hist.belief_std_Cthr = zeros(maxIter, 4, scales);
        hist.Cdiv = zeros(maxIter, 4, scales);
        hist.RMSEyest = zeros(1, maxIter);
        hist.x_mse = zeros(1, maxIter);
        hist.r1 = zeros(nx,ny,maxIter);
        hist.C_thr = zeros(nx,ny,maxIter);
    end
    
    
    
    C_damp= 0; % for damping
    time_init = tic;   
    
    for iter = 1:maxIter
        
        if saveHist
            hist.r1(:,:,iter) = r;
        end
        
        if saveHist
            hist.RMSEyest(iter) = sqrt(mean(y_cov(:)));
        end
        
        if verbose % band-wise mse of unthresholded image
            disp(['                 ITERATION ', num2str(iter)])
            fprintf( ...
                'K-space message: true RMSE = %f\n', ...
                sqrt(mean(abs(r(:)-x0(:)).^2)));
            for s=1:scales
                subbands = fieldnames(C{s});
                fprintf(['\tAt scale ' int2str(s) '\n']);
                for i = 1:numel(subbands)
                    band_name = subbands{i};
                    fprintf( ...
                        '\t\tIn band %s: true RMSE = %f\tmessage std = %f\n', ...
                        band_name, ...
                        sqrt(mean(abs(C{s}.(band_name)(:)-W0{s}.(band_name)(:)).^2)), ...
                        sqrt(tau{s}.(band_name)));
                end
            end
        end
        
        % thresholding
        if sureFlag
            [C_thr_new, err, df] = multiscaleSUREsoft(C, tau);
        else      
            [C_thr_new, err, df] = multiscaleComplexSoft(C, tau, lambda);
        end
            
        
        if saveHist
            hist.C_thr(:,:, iter) = pyramid(C_thr_new);
        end
             
        C_thr = C_thr_new;
        
        % calculate x estimate and MSE via unweighted density compensatoin
        xk_tilde = multiscaleRecon(C_thr);
        x_tilde_ft = fftnc(xk_tilde);
        z2 = mask.*(dcoil - x_tilde_ft);
        xk = ifftnc(x_tilde_ft + z2);
        hist.x_mse(iter) = immse(xk,x0);
        
        hist.timer(iter) = toc(time_init);
        if hist.timer(iter) > maxTime
            break
        end
        
        %Onsager correction in wavelet domain the OAMP/VAMP way
        for s=1:scales
            subbands = fieldnames(C_thr{s});           
            for i = 1:numel(subbands)
                band_name = subbands{i};
                
                alpha{s}.(band_name) = mean(df{s}.(band_name)(:))/2;   
                C_tilde{s}.(band_name) =  (C_thr{s}.(band_name) - alpha{s}.(band_name)*C{s}.(band_name));
                
                if denoiserDiv == 1    
                    % VDAMP-S
                    Cdiv = C{s}.(band_name)(:)'*(C_tilde{s}.(band_name)(:))./norm(C_tilde{s}.(band_name)(:))^2;
                else
                    % VDAMP-alpha
                    Cdiv = 1./(1 - alpha{s}.(band_name));
                end
                
                C_tilde{s}.(band_name) = C_tilde{s}.(band_name).*Cdiv;

                if saveHist         
                    hist.true_err_C(iter,i,s) = sqrt(mean(abs(C{s}.(band_name)(:)-W0{s}.(band_name)(:)).^2));
                    hist.belief_std_C(iter,i,s) = sqrt(tau{s}.(band_name));
                    hist.true_err_Cthr(iter,i,s) = sqrt(mean(abs(C_thr{s}.(band_name)(:)-W0{s}.(band_name)(:)).^2));
                    hist.belief_std_Cthr(iter,i,s) = sqrt(mean(abs(C_thr{s}.(band_name)(:)-C{s}.(band_name)(:)).^2)-tau{s}.(band_name) + 2*err{s}.(band_name));
                    hist.Cdiv(iter,i,s) = Cdiv;
                end        
            end
        end
 
        if verbose % band-wise progress of thresholded image 
             x = multiscaleRecon(C_thr);
            fprintf('Wavelet belief: true RMSE = %f\n', sqrt(mean(abs(x(:)-x0(:)).^2)));
            for s=1:scales
                subbands = fieldnames(C{s});
                fprintf(['\tAt scale ' int2str(s) '\n']);
                for i = 1:numel(subbands)
                    band_name = subbands{i};
                    fprintf( ...
                        '\t\tIn band %s: true RMSE = %f\tbelief std = %f\tSURE = %f\n', ...
                        band_name, ...
                        sqrt(mean(abs(C_thr{s}.(band_name)(:)-W0{s}.(band_name)(:)).^2)), sqrt(err{s}.(band_name)),...
                        sqrt(mean(abs(C_thr{s}.(band_name)(:)-C{s}.(band_name)(:)).^2)-tau{s}.(band_name) + 2*err{s}.(band_name)));
                end
            end
        end
      
        % ****Reweighted gradient step****
        r_tilde = multiscaleRecon(C_tilde);
        z = mask.*(dcoil - fftnc(r_tilde));
        
        r = r_tilde + ifftnc(inv_p .* z);    
        C = multiscaleDecomp(r, scales);
                
        % coloured noise power spectrum re-estimation
        y_cov = mask.* (inv_p_m1 .* inv_p .* abs(z).^2+ inv_p.*var0);
        for s = 1:scales % to wavelet domain
            tau{s}.H = specX(:, s, 2)'* y_cov * specY(:, s, 1); 
            tau{s}.V = specX(:, s, 1)'* y_cov * specY(:, s, 2); 
            tau{s}.D = specX(:, s, 2)'* y_cov * specY(:, s, 2);
        end
        tau{s}.A =  specX(:, s, 1)'* y_cov * specY(:, s, 1); 
    end
  
    if ~saveHist
        hist.r1 = r;
    end
 
    x_hat = xk;
end


