% generates and saves sista_weights using the power iteration method
% for various masks and dimensions
%
% The Software does not have 510K approval,  U.S. Food and Drug Administration (FDA) clearance,
% European Medicines Agency (EMA) clearance or CE mark. The Software is intended research purposes only,
% not for diagnostic or clinical use.
%
% Copyright (C) 2019  Charles Millard
% Copyright (C) 2019  Siemens Healthineers

clear all; 
addpath(genpath(pwd))
%%

rng('default');

%%
scales = 4;

maxIter = 100; % run power iteration mthod for 100 iterations

for target_delta = 1./[4,6,8]
    rng(811);
    nx = 256; ny = 256;
    prob_map = genPDF([nx,ny], 8, target_delta, 2, 0,0);
    mask = binornd(1, prob_map, nx,ny);
    sista_div = sista_div_generator(scales, mask, maxIter);
    save(['/home/user/VDAMP-OJSP/sista_weights/256x256/us_fac', num2str(round(1/target_delta))], 'sista_div');
    disp([num2str([nx,ny, round(1./target_delta)]), 'x complete'])
    

    rng(811);
    nx = 208; ny = 416; % for cardiac MR image
    prob_map = genPDF([nx,ny], 8, target_delta, 2, 0,0);
    mask = binornd(1, prob_map, nx,ny);
    sista_div = sista_div_generator(scales, mask, maxIter);
    
    save(['/home/user/VDAMP-OJSP/sista_weights/208x416/us_fac', num2str(round(1/target_delta))], 'sista_div');
    disp([num2str([nx,ny, round(1./target_delta)]), 'x complete'])
end

for target_delta = 1./[4,6,8,10,12] % large images (plus shepp-logan
    rng(811);
    nx = 512; ny = 512;
    prob_map = genPDF([nx,ny], 8, target_delta, 2, 0,0);
    mask = binornd(1, prob_map, nx,ny);
    sista_div = sista_div_generator(scales, mask, maxIter);
    
    save(['/home/user/VDAMP-OJSP/sista_weights/512x512/us_fac', num2str(round(1/target_delta))], 'sista_div');
    disp([num2str([nx,ny, round(1./target_delta)]), 'x complete'])
end

function sista_div = sista_div_generator(scales, mask, maxIter)
%generates the SISTA divisors using the power iteration method for Haar
%wavelets for a given mask

for s =1:scales
    sista_div{s} = struct('H', 0, 'V', 0, 'D', 0);
end

sista_div{scales}.A = 0;

subb_count = 0;

for s2= 1:scales
    subbands2 = fieldnames(sista_div{s2});
    for i2 = 1:numel(subbands2)
        band_name2 = subbands2{i2};
        rho = zeros(scales,4); % maximum eigenvalues
        for s1 =1:scales
            subbands1 = fieldnames(sista_div{s1});           
            for i = 1:numel(subbands1)
                band_name1 = subbands1{i};

                % random initialization
                b = multiscaleDecomp(zeros(size(mask)), scales);
                b{s1}.(band_name1) = normrnd(0,1, size(b{s1}.(band_name1)));
                
                az = multiscaleDecomp(zeros(size(mask)), scales);
                qz = multiscaleDecomp(zeros(size(mask)), scales);
                for iter = 1:maxIter                
                    a =  multiscaleDecomp(ifftnc(mask.*fftnc(multiscaleRecon(b))), scales);           
                    az{s2}.(band_name2) = a{s2}.(band_name2);

                    q = multiscaleDecomp(ifftnc(mask.*fftnc(multiscaleRecon(az))), scales);                 
                    qz{s1}.(band_name1) = q{s1}.(band_name1);
             
                    if iter == maxIter
                        z_py_last = pyramid(qz);
                    else
                        z_py = pyramid(qz);
                        qz_normed = z_py./norm(z_py(:));
                        b = pyramidInv(qz_normed, scales);
                    end
                end 
                rho(s1,i) = qz_normed(:)'*z_py_last(:)./(qz_normed(:)'*qz_normed(:));
            end
        end
    sista_div{s2}.(band_name2) = sum(sqrt(real(rho(:))));
    subb_count = subb_count + 1;
    end
end

end
