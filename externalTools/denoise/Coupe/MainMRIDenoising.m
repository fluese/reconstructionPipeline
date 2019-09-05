% Pierrick Coupe - pierrick.coupe@gmail.com
% Jose V. Manjon - jmanjon@fis.upv.es
% LaBRI UMR 5800
% Universite de Bordeaux 1
%
% Copyright (C) 2008-2013 Pierrick Coupe and Jose V. Manjon
% ***************************************************************************/

warning off;
clc;
clf;
clear all;
close all

%install

fprintf('%s \n', 'Welcome to MRI Denoising software:')
fprintf('%s \n\n', 'Copyright (C) 2008-2013 Pierrick Coupe and Jose Manjon')

verbose =0;

addpath 'spm8';
addpath 'gui';
packagepath = pwd;
addpath(genpath(fullfile(packagepath, 'MRIDenoisingPackage')));


[namein, pathin, filterindex] = uigetfile({  '*.nii','NIFTI image (*.nii)'}, 'Select one or several files (using +CTRL or +SHIFT)','MultiSelect', 'on');
if isequal(namein,0) | isequal(pathin,0)
    disp('User pressed cancel')
else
    
    
    [flag method beta patchradius searchradius suffixfile rician] = gui;
    flag
    method
    beta
    patchradius
    searchradius
    suffixfile
    rician
    
    fprintf('%s\n\n', 'Selected parameters')
    if (method==1)
        com = sprintf('Optimized NLM Filter\n');
        disp(com)
    end
    if (method==2)
        com = sprintf('Adaptative ONLM Filter\n');
        disp(com)
    end
    if (method==3)
        com = sprintf('Multi-resolution ONLM Filter\n');
        disp(com)
    end
    
    if (method==4)
        com = sprintf('Oracle DCT Filter\n');
        disp(com)
    end
    
    if (method==5)
        com = sprintf('Prefiltered rationally invariant NLM Filter\n');
        disp(com)
    end
    
    com = sprintf('beta: %0.1f', beta);
    disp(com)
    com = sprintf('patch size: %dx%dx%d voxels', 2*patchradius+1, 2*patchradius+1, 2*patchradius+1);
    disp(com)
    com = sprintf('search volume size: %dx%dx%d voxels', 2*searchradius+1, 2*searchradius+1, 2*searchradius+1);
    disp(com)
    if (rician==1)
        com = sprintf('Rician noise model\n');
        disp(com)
    else
        com = sprintf('Gaussian noise model\n');
        disp(com)
    end
    
    
    if (flag==1)
        
        if (iscell(namein))
            nbfile = size(namein,2);
        else
            nbfile = 1;
        end
        
        for f = 1 : nbfile
            
            
            
            if(nbfile>1)
                filenamein = namein{f};
            else
                filenamein = namein;
            end
            
            disp(['Input file : ', fullfile(pathin, filenamein)])
            [pathstr, name_s, ext]=fileparts(fullfile(pathin, filenamein));
            nout=[name_s suffixfile ext];
            pathout = pathin;
            
            disp(['Output file: ', fullfile(pathout, nout)])
            
            
            
            cname = fullfile(pathin, filenamein);
            VI=spm_vol(cname);
            ima=spm_read_vols(VI);
            s=size(ima);
            
            % fixed range
            map = isnan(ima(:));
            ima(map) = 0;
            map = isinf(ima(:));
            ima(map) = 0;
            mini = min(ima(:));
            ima = ima - mini;
            maxi=max(ima(:));
            ima=ima*256/maxi;
            
            
            %Noise estimation
            if(method~=2)
                verbose =0;
                [hfinal, ho, SNRo, hbg, SNRbg] = MRINoiseEstimation(ima, rician, verbose);
                
                if(isnan(hfinal))
                    disp('error during noise estimation')
                    continue;
                end
            end
            
            
            verbose =1;
            
            % Denoising
            if(method==1)
                MRIdenoised = MRIDenoisingONLM(ima, hfinal, beta, patchradius,  searchradius, rician, verbose);
            end
            
            if(method==2)
                MRIdenoised = MRIDenoisingAONLM(ima, patchradius,  searchradius, beta, rician, verbose);
            end
            
            if(method==3)
                MRIdenoised = MRIDenoisingMRONLM(ima, hfinal, beta, patchradius,  searchradius, rician, verbose);
            end
            
            if(method==4)
                MRIdenoised = MRIDenoisingODCT(ima, hfinal, beta, rician, verbose);
            end
            
            if(method==5)
                MRIdenoised = MRIDenoisingPRINLM(ima, hfinal, beta, rician, verbose);
            end
            map = find(MRIdenoised<0);
            MRIdenoised(map)=0;
            
            
            
            % Original intensity range
            MRIdenoised= MRIdenoised*maxi/256;
            MRIdenoised =MRIdenoised + mini;
            
            
            VO = VI; % copy input info for output image
            
            outfilename =fullfile(pathout, nout);
            VO.fname = outfilename;
            VO.dim=s;
            spm_write_vol(VO,MRIdenoised(:,:,:));
            
        end
    end
end


