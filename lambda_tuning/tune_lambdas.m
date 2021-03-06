% Tune sparse weighting of FISTA and SISTA with an exhaustive search
%
%
% The Software does not have 510K approval,  U.S. Food and Drug Administration (FDA) clearance,
% European Medicines Agency (EMA) clearance or CE mark. The Software is intended for research purposes only,
% not for diagnostic or clinical use.
%
% Copyright (C) 2020  Charles Millard
% Copyright (C) 2020  Siemens Healthineers

clear all;

all_ims = {'shepplogan', 'brain', 'cardiac', 'barbara', 'boat',  'cameraman', 'house', 'peppers'};
N_im = numel(all_ims);

for im_idx = 1:N_im
    im_type = all_ims{im_idx};
  
    disp(all_ims{im_idx})
    if strcmp(all_ims{im_idx}, 'shepplogan')
        us_fac = 1./[8,10,12];
    else
        us_fac = 1./[4,6,8];
    end

    for ii = 1:numel(us_fac)
        [fista_opt, sista_opt] = exh_search(im_type, us_fac(ii));
 
        fista_opt_all.(im_type).(['us', num2str(round(1/us_fac(ii)))]) = fista_opt;
        sista_opt_all.(im_type).(['us', num2str(round(1/us_fac(ii)))]) = sista_opt;
    end  
end

save('fista_opt', 'fista_opt_all');
save('sista_opt', 'sista_opt_all');

function [fista_opt, sista_opt] = exh_search(im_type, target_delta)
    % finds best lambda for im_type at target_delta = n/N
    rng('default');
    rng(811);

    switch im_type
        case 'shepplogan'
            x0 = phantom('modified shepp-logan', 512);
        otherwise
            x0 = imread([num2str(im_type), '.png']);
    end

    x0 = double(rescale(x0));
    [nx, ny] = size(x0);    
    prob_map = genPDF([nx,ny], 8, target_delta, 2, 0,0);

    %% generate data
    SNR = 40;

    mask = binornd(1, prob_map, nx,ny);

    var0 = mean(abs(x0(:)).^2)/(10^(SNR/10));    
    noise = normrnd(0, sqrt(var0), nx,ny)./sqrt(2) + 1i* normrnd(0, sqrt(var0), nx,ny)./sqrt(2);

    dcoil = mask.*(fftnc(x0) + noise);

    %% set algorithm options

    opts.verbose =0;
    opts.fistaFlag = 2;
    opts.rho = 1;
    opts.saveHist = 0;
    opts.scales = 4;
    opts.maxTime = 10000;

    load(['sista_weights/', num2str(nx), 'x', num2str(ny), '/us_fac', num2str(round(1/target_delta))]);
    %% optional exhaustive search for optimal lambda 


    opts.SURE = 0;
    opts.maxIter = 100;
    trial_lambda = [1:40]; % exhaustive search range
    disp('Tuning FISTA/SISTA...')
    x_mse_trial_fista = zeros(1,numel(trial_lambda));
    x_mse_trial_sista = zeros(1,numel(trial_lambda));
    for l = 1:numel(trial_lambda)
        opts.lambda = trial_lambda(l);
        opts.sista_div = 1;
        x_mse_trial_fista(l) = immse(FISTA(dcoil, mask, x0, opts),x0);
        opts.sista_div = sista_div;
        x_mse_trial_sista(l) = immse(FISTA(dcoil, mask, x0, opts),x0);
    end
    opt_idx_fista = find(x_mse_trial_fista == min(x_mse_trial_fista));
    opt_idx_sista = find(x_mse_trial_sista == min(x_mse_trial_sista));

    disp(['FISTA lambda: ',num2str(trial_lambda(opt_idx_fista))]);
    disp(['SISTA lambda: ',num2str(trial_lambda(opt_idx_sista))]);

    fista_opt = trial_lambda(opt_idx_fista);
    sista_opt = trial_lambda(opt_idx_sista);
end
