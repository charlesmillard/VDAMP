% Comparion between VDAMP-alpha, VDAMP-S, FISTA, S-FISTA and SURE-IT for uniform,
% two-level or polynomial undersampling. 
%
% Uses pre-tuned sparse weighting lambdas in lambda_tuning folder. Also
% uses pre-determined per-subband sista weights from sista_weights folder.
%
% Note that to run any reconstruction tasks, it is necessary to use
% an image type and acceleration factor N/n used in the paper: 4,6,8 for
% all images except the shepp-logan, which uses 8, 10 12. If you would like
% to run FISTA/SISTA with other masks, you can use the scipts in lambda_tuning 
% and sista_weights to select them
%
% The Software does not have 510K approval,  U.S. Food and Drug Administration (FDA) clearance,
% European Medicines Agency (EMA) clearance or CE mark. The Software is intended for research purposes only,
% not for diagnostic or clinical use.
%
% Copyright (C) 2020  Charles Millard
% Copyright (C) 2020  Siemens Healthineers

clear all; 
addpath(genpath(pwd))
%%

rng('default');
rng(811);

%% im_type options are 'shepplogan' or those listed in test_images folder

im_type = 'cardiac';
target_delta = 1/4; % n/N

%%
switch im_type
    case 'shepplogan'
        x0 = phantom('modified shepp-logan', 512);
    otherwise
        x0 = imread([num2str(im_type), '.png']);
end

x0 = double(rescale(x0));
[nx, ny] = size(x0);

%% lustig's VDS scheme 

prob_map = genPDF([nx,ny], 8, target_delta, 2, 0,0);


%% generate data
SNR = 40;
        
mask = binornd(1, prob_map, nx,ny);
delta=mean(mask(:));

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

%% load pre-tuned sista weights
try
    load(['sista_weights/', num2str(nx), 'x', num2str(ny), '/us_fac', num2str(round(1/target_delta))]);
catch
    error('SISTA weights not yet tuned of this undersampling factor')
end
    
opts.sista_div = sista_div;

%% load pre-tuned lambda
load(['lambda_tuning/fista_opt'])
load(['lambda_tuning/sista_opt'])
    
try
    lambda_ista = fista_opt_all.(im_type).(['us', num2str(round(1/target_delta))]);
catch
    error('this fista lambda has not been pre-tuned')
end

try
    lambda_sista = sista_opt_all.(im_type).(['us', num2str(round(1/target_delta))]);
catch
    error('this sista lambda has not been pre-tuned')
end

%% Run FISTA and VDAMP
if strcmp(im_type, 'shepplogan')
    opts.maxIter = 1000;
else 
    opts.maxIter = 100;
end

% FISTA
opts.lambda = lambda_ista;    
opts.SURE = 0;
opts.sista_div = 1;
disp('Running FISTA...')
[x_ista, hist_ista]  = FISTA(dcoil, mask, x0, opts);

% S-FISTA
opts.lambda = lambda_sista;
opts.SURE = 0;
opts.sista_div = sista_div;
disp('Running SISTA...')
[x_sista, hist_sista]  = FISTA(dcoil,  mask, x0, opts);

% SURE-IT
opts.SURE = 1;
opts.sista_div = 1;
disp('Running SURE-IT...')
[x_sureit, hist_sureit] = FISTA(dcoil,  mask, x0, opts);

% VDAMP-alpha
opts.denoiserDiv = 2;
disp('Running VDAMP...')
[x_vdamp,hist_vdamp]  = VDAMP(dcoil, mask, prob_map, var0, x0, opts);

% VDAMP-S
opts.denoiserDiv = 1;
disp('Running VDAMP-S...')
[x_vdamp_s,hist_vdamp_s]  = VDAMP(dcoil, mask, prob_map, var0, x0, opts);
%%
ista_res = abs(x0 - x_ista);
sista_res = abs(x0 - x_sista);
sureit_res = abs(x0 - x_sureit);
vdamp_res = abs(x0 - x_vdamp);
vdamp_s_res = abs(x0 - x_vdamp_s);

thr = 0.3*max([sureit_res(:); ista_res(:)]); %clip for visualization

ista_res(ista_res>thr) = thr;
sista_res(sista_res>thr) = thr;
sureit_res(sureit_res>thr) = thr;
vdamp_res(vdamp_res>thr) = thr;
vdamp_s_res(vdamp_s_res>thr) = thr;


figure('Name', 'FISTA and VDAMP'); 
subplot(2,6,1); 
imagesc(x0);
axis image off;
colormap(gca,'gray')
set(gca,'FontName','times')
title('x0')

subplot(2,6,2); 
imagesc(abs(x_ista), [0 1]);
axis image off;
colormap(gca,'gray')
set(gca,'FontName','times')
title('FISTA')

subplot(2,6,3); 
imagesc(abs(x_sista), [0 1]);
axis image off;
colormap(gca,'gray')
set(gca,'FontName','times')
title('SISTA')

subplot(2,6,4); 
imagesc(abs(x_sureit), [ 0 1]);
axis image off;
colormap(gca,'gray')
set(gca,'FontName','times')
title('SURE-IT')

subplot(2,6,5); 
imagesc(abs(x_vdamp), [ 0 1] );
axis image off;
colormap(gca,'gray')
set(gca,'FontName','times')
title('VDAMP-\alpha')

subplot(2,6,6); 
imagesc(abs(x_vdamp_s), [ 0 1] );
axis image off;
colormap(gca,'gray')
set(gca,'FontName','times')
title('VDAMP-S')

subplot(2,6,7);
imagesc(mask);
axis image off;
colormap(gca,'gray')
set(gca,'FontName','times')
title('Mask')

subplot(2,6,8);
imagesc(ista_res, [0 thr]);
axis image off;
colormap(gca,'jet')
set(gca,'FontName','times')
title('FISTA error')

subplot(2,6,9);
imagesc(sista_res, [0 thr]);
axis image off;
colormap(gca,'jet')
set(gca,'FontName','times')
title('FISTA error')

subplot(2,6,10);
imagesc(sureit_res, [0 thr]);
axis image off;
colormap(gca,'jet')
set(gca,'FontName','times')
title('SURE-IT error')

subplot(2,6,11);
imagesc(vdamp_res, [0, thr ])
axis image off;
colormap(gca,'jet')
set(gca,'FontName','times')
title('VDAMP error')

subplot(2,6,12);
imagesc(vdamp_s_res, [0, thr ])
axis image off;
colormap(gca,'jet')
set(gca,'FontName','times')
title('VDAMP S error')

%%
x_mnsq = mean(abs(x0(:)).^2);
NMSE = @(x_mse) 10*log10(x_mse./x_mnsq);

line_colour =  [0.1, 0.6, 0; 1, 0.8, 0; 1, 0, 0; 0, 0.5, 1; 0,0,0];

figure('Renderer', 'painters', 'Position', [10 10 200 200]);
hold on;
plot([0: numel(hist_ista.x_mse)-1 ], [NMSE(hist_ista.x_mse)], '--', 'LineWidth', 1.5, 'Color', line_colour(1,:));
plot([0: numel(hist_sista.x_mse)-1 ], [NMSE(hist_sista.x_mse)], '--', 'LineWidth', 1.5, 'Color', line_colour(5,:));
plot([0: numel(hist_sureit.x_mse)-1 ], [NMSE(hist_sureit.x_mse)], '--', 'LineWidth', 1.5, 'Color',line_colour(2,:));
plot([0: numel(hist_vdamp.x_mse)-1 ], [NMSE(hist_vdamp.x_mse)], '-', 'LineWidth', 1.5, 'Color',line_colour(3,:));
plot([0: numel(hist_vdamp_s.x_mse)-1 ], [NMSE(hist_vdamp_s.x_mse)], '-', 'LineWidth', 1.5, 'Color',line_colour(4,:));

xlabel(['Iteration k']);

ylabel('NMSE (dB)');

set(gca, 'FontName', 'Times' );
grid on

legend('FISTA', 'S-FISTA', 'SURE-IT','VDAMP-\alpha', 'VDAMP-S');
title(['N/n =', num2str(round(1/delta))])
hold off
%%

figure('Renderer', 'painters', 'Position', [10 10 200 200]);
hold on;
plot(hist_ista.timer, [NMSE(hist_ista.x_mse)], '--', 'LineWidth', 1.5, 'Color', line_colour(1,:));
plot(hist_sista.timer, [NMSE(hist_sista.x_mse)], '--', 'LineWidth', 1.5, 'Color', line_colour(5,:));
plot(hist_sureit.timer, [NMSE(hist_sureit.x_mse)], '--', 'LineWidth', 1.5, 'Color',line_colour(2,:));
plot(hist_vdamp.timer, [NMSE(hist_vdamp.x_mse)], '-', 'LineWidth', 1.5, 'Color',line_colour(3,:));
plot(hist_vdamp_s.timer, [NMSE(hist_vdamp_s.x_mse)], '-', 'LineWidth', 1.5, 'Color',line_colour(4,:));

xlabel(['Time (s)']);

xlim([0,  5]);
ylabel('Time (s)');
set(gca, 'FontName', 'Times' );
grid on

legend('FISTA', 'S-FISTA', 'SURE-IT','VDAMP-\alpha', 'VDAMP-S');
title(['N/n =', num2str(round(1/delta))])
hold off

%% iters until within 1dB of convergent NMSE

conv_thr = 0.1;

conv_ista = find(abs(NMSE(hist_ista.x_mse) - NMSE(hist_ista.x_mse(end)))<conv_thr ,1);
conv_sista = find(abs(NMSE(hist_sista.x_mse) - NMSE(hist_sista.x_mse(end)))<conv_thr ,1);
conv_sureit = find(abs(NMSE(hist_sureit.x_mse) - NMSE(hist_sureit.x_mse(end)))<conv_thr ,1);
conv_vdamp = find(abs(NMSE(hist_vdamp.x_mse) - NMSE(hist_vdamp.x_mse(end)))<conv_thr ,1);
conv_vdamp_s = find(abs(NMSE(hist_vdamp_s.x_mse) - NMSE(hist_vdamp_s.x_mse(end)))<conv_thr ,1);

w0 = multiscaleDecomp(x0, opts.scales);

disp('** FISTA **')
[re, im] = subband_kurtosis(w0, hist_ista.r1(:,:,end));
disp(['NMSE: ', num2str(NMSE(hist_ista.x_mse(end)))])
disp(['Convergent index/time: ', num2str([conv_ista, hist_ista.timer( conv_ista)])])
disp(['Kurtosis re/im: ' , num2str([re, im])])

disp('** SISTA **')
[re, im] =  subband_kurtosis(w0, hist_sista.r1(:,:,end));
disp(['NMSE: ', num2str(NMSE(hist_sista.x_mse(end)))])
disp(['Convergent index/time: ', num2str([conv_sista, hist_sista.timer( conv_sista)])])
disp(['Kurtosis re/im: ' , num2str([re, im])])

disp('** SUREIT **')
[re, im] = subband_kurtosis(w0, hist_sureit.r1(:,:,end));
disp(['NMSE: ', num2str(NMSE(hist_sureit.x_mse(end)))])
disp(['Convergent index/time: ', num2str([conv_sureit, hist_sureit.timer( conv_sureit)])])
disp(['Kurtosis re/im: ' , num2str([re, im])])

disp('** VDAMP **')
[re, im] = subband_kurtosis(w0, hist_vdamp.r1(:,:,end));
disp(['NMSE: ', num2str(NMSE(hist_vdamp.x_mse(end)))])
disp(['Convergent index/time: ', num2str([conv_vdamp, hist_vdamp.timer( conv_vdamp)])])
disp(['Kurtosis re/im: ' , num2str([re, im])])

disp('** VDAMP-S **')
[re, im]  = subband_kurtosis(w0, hist_vdamp_s.r1(:,:,end));
disp(['NMSE: ', num2str(NMSE(hist_vdamp_s.x_mse(end)))])
disp(['Convergent index/time: ', num2str([conv_vdamp_s, hist_vdamp_s.timer( conv_vdamp_s)])])
disp(['Kurtosis re/im: ' , num2str([re, im])])


