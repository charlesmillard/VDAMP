% Detailed evaluation of only VDAMP. Feel free to change the image type,
% undersampling factor and algorithm opts
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

im_type = 'brain'; 

target_delta = 1/5; % n/N

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

%% VDAMP options
opts.maxIter = 30;
opts.maxTime = 100;
opts.verbose =0 ;
opts.scales = 4;
opts.SURE = 1;
opts.saveHist = 1;
opts.wavType = 'db1';
opts.denoiserDiv = 1; % 1 for VDAMP-S or 2 for VDAMP-alpha

%% Run VDAMP
t = tic;
[x_hat, hist]  = VDAMP(dcoil,  mask, prob_map, var0, x0, opts);
timer = toc(t);


%% How well is the error modelled by tau? 

w0 = multiscaleDecomp(x0, opts.scales, opts.wavType);

% ground truth subband norms for later NMSE calculation
for s = 1:opts.scales
    subbands = fieldnames(w0{s});
    for i = 1:numel(subbands)
        band_name = subbands{i};
        C0norm(i,s) = norm(w0{s}.(band_name));
        C0sz(i,s) = numel(w0{s}.(band_name));
    end
end

k_idx = [1:opts.maxIter]-1; 
col_arr = [0,0.5,1; 1, 0, 0; 0.5, 0, 0.5; 0, 0.5 0]; %line colours

figure('Name', 'Tracking the subband-wise NMSE: lines are ground truth and crosses are the estimate'); 
for s = 1:opts.scales
    %subplot(2,ceil(opts.scales/2),s); hold on;
    subplot(1,ceil(opts.scales),s); hold on;
    
    title(['Scale ', num2str(s)]);
    set(gca, 'FontName','times')
    
    ylabel('NMSE (dB)');
    xlabel('Iteration');    
 
    if s == opts.scales
        detail_max = 5;
    else
        detail_max = 4;
    end
    
    detail=1;
    col = 1;
    while detail<detail_max
        line_colour = col_arr(col,:);
        plot(k_idx, 10*log10(C0sz(detail,s)*squeeze(hist.true_err_C(:,detail,s)).^2./C0norm(detail,s).^2), '-', 'Color' , line_colour)
        plot(k_idx, 10*log10(C0sz(detail,s)*squeeze(hist.belief_std_C(:,detail,s)).^2./C0norm(detail,s).^2),'x', 'HandleVisibility','on', 'Color', line_colour);
        detail=detail+1;
        col = col+1;
    end
    
    hold off;
end

legend({'True Horizontal', 'Horizontal Model', 'True Vertical', 'Vertical Model',  'True Diagonal', 'Diagonal Model', 'True Coarse', 'Coarse Model'}, 'NumColumns', 4)
%% Visualise difference in wavelet domain. 

I0 = pyramid(w0);
it_choice = [2, 3, 4]; % which iterations to look at

%
figure('Name', 'Difference between r and ground truth');
for iter = 1:numel(it_choice)  
    subplot(1,numel(it_choice),iter);
    ii = it_choice(iter);
    
    C = multiscaleDecomp(hist.r1(:,:,ii),opts.scales, opts.wavType);
    IC = pyramid(C);
    diff = abs(IC-I0);
    
    if iter == 1
        thr = 0.3; % choose clip
        diff(diff>thr)= thr;
        diff = rescale(diff);
    else
        diff(diff>thr) = thr;
        diff = rescale(diff)*max(diff(:))/thr;
    end

    
    imagesc(diff, [0, 1]); colormap jet;
    title(['k=', num2str(ii-1)]);
    colorbar off;
    axis image off;
    set(gca,'FontName','times')
end

%% QQ plots - is the effective noise Gaussian? 

it_choice = [1, 6, 21]; % which iterations to look at

div = 1; %normalisation
qqylim = [-1, 1];


figure('Name', 'QQ plots');
for iter = 1:numel(it_choice)
    subplot(3,numel(it_choice),iter);
    ii = it_choice(iter);
    C = multiscaleDecomp(hist.r1(:,:,ii),opts.scales, opts.wavType);
    
    diff = (w0{1}.D - C{1}.D)/div; 
    
    p = qqplot(real(diff(:)));
    set(p, 'Marker', '.');
    
    xh = xlabel('Standard normal quantiles');
    set(xh, 'FontName', 'times');
    set(gca,'FontSize',8)
    
    yh = ylabel('Quantiles of input sample');
    set(yh, 'FontName', 'times');
    
    axis square;
    title(['Diagonal, scale 1, k=', num2str(ii-1)], 'FontWeight','bold', 'FontName', 'Times new roman');
    ylim([-1 1]);
    xlim([-5, 5]);
end

for iter = 1:numel(it_choice) 
    subplot(3,numel(it_choice),iter+numel(it_choice));
    ii = it_choice(iter);
    C = multiscaleDecomp(hist.r1(:,:,ii),opts.scales, opts.wavType);
    
    diff = (w0{2}.H - C{2}.H)/div; 
    
    p = qqplot(imag(diff(:)));
    set(p, 'Marker', '.');
    
    xh = xlabel('Standard normal quantiles');
    set(xh, 'FontName', 'times');
    
    yh = ylabel('Quantiles of input sample');
    set(yh, 'FontName', 'times');
    axis square;
    set(gca,'FontSize',8)
    
    title(['Horizonal, scale 2, k=', num2str(ii-1)], 'FontWeight','bold', 'FontName', 'Times new roman');
    ylim([-1 1]);
    xlim([-5, 5]);
end

for iter = 1:numel(it_choice)
    subplot(3,numel(it_choice),iter+2*numel(it_choice));
    ii = it_choice(iter);
    C = multiscaleDecomp(hist.r1(:,:,ii),opts.scales, opts.wavType);
    
    diff = (w0{4}.V - C{4}.V)/div; 
    
    p = qqplot(real(diff(:)));
    set(p, 'Marker', '.');
    xh = xlabel('Standard normal quantiles');
    set(xh, 'FontName', 'times');
    
    yh = ylabel('Quantiles of input sample');
    set(yh, 'FontName', 'times');
    axis square;
    set(gca,'FontSize',8)
    
    title(['Vertical, scale 4, k=', num2str(ii-1)], 'FontWeight','bold', 'FontName', 'Times new roman');
    ylim([-1 1]);
    xlim([-5, 5]);
end


%% view reconstructed image

figure('Name', 'VDAMP reconstruction and error'); 
subplot(2,2,1);
imagesc(x0); 
colormap(gca, 'gray')
axis image off;
title('x0');

subplot(2,2,2)
imagesc(abs(x_hat), [0 max(x0(:))]);
colormap(gca, 'gray')
axis image off;
title('VDAMP')

subplot(2,2,3)
imshow(abs(ifftnc(dcoil)), [0 max(x0(:))]);
title('Zero-filled')

subplot(2,2,4)
imagesc(abs(x_hat-x0)); 
axis image off;
colormap(gca, 'jet')
colorbar off;
title('VDAMP error')

%%

disp(['Undersampling factor: ', num2str(1/delta)])

x_mnsq = mean(abs(x0(:)).^2);
NMSE = @(x_mse) 10*log10(x_mse./x_mnsq);

disp(['NMSE is: ', num2str(NMSE(immse(x0,x_hat)))])

%% Subband-wise kurtosis

[re, im] = subband_kurtosis(w0, hist.r1(:,:,end), opts.wavType);
disp(['Kurtosis re/im: ' , num2str([re, im])])