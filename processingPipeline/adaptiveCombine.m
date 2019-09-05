function [recon,wfull,norm]=adaptiveCombine(settings, im)
%   Adaptive recon based on Walsh et al.
%   Walsh DO, Gmitro AF, Marcellin MW.
%   Adaptive reconstruction of phased array MR imagery.
%   Magn Reson Med. 2000 May;43(5):682-90.
%
%    and
%
%   Mark Griswold, David Walsh, Robin Heidemann, Axel Haase, Peter Jakob.
%   The Use of an Adaptive Reconstruction for Array Coil Sensitivity Mapping and Intensity Normalization,
%   Proceedings of the Tenth  Scientific Meeting of the International Society for Magnetic Resonance in Medicine pg 2410 (2002)
%
%   Changes by Falk Luesebrink (14.12.2017)

%%
% Setting up variables.
cnt     = 0;
sz = size(im);
n  = ones(1,3);
n(1:numel(sz)-1) = sz(2:end);

%%
% Enables singular value decomposition, referred to as coil compression. It
% is accelerates the coil combination process as less coils have to be
% considered and is supposed to be much better for phase correction.
%
% Disabled by default.
if settings.options.xspace.coilCompression
    nc_svd = min(max(16,floor(settings.parameters.nCh/2)),settings.parameters.nCh);
else
    nc_svd = settings.parameters.nCh;
end

%%
% Enables processing slice by slice. Using 3D data may produce better
% background suppression and, therefore, higher SNR. See paragraph for
% setting the blocksize for in depth explation.
%
% Disabled by default.

%%
% This section allows for intensity normalization suggested in the
% publication. From what I have seen the center of the image seems to be
% intensified, while the rest is suppressed. I haven't seen any image that
% looked overall better than without this normalization.
%
% There may, however, be better approaches to normalize the image than the
% one used here.
%
% Disabled by default.
%
if ~isfield(settings.parameters,'donorm')
    settings.parameters.donorm = false;
end

%%
% The original publication suggests to use a noise correlation matrix for
% better performance.
%
% In the current implementation the coils are decorrelated before running
% this script by making use of a noise reference scan. See according
% section in the OTHER scripts (insert name once everything is finalized).
%
% As it is de-correlated already, we can simply use the unity matrix.
%
if ~isfield(settings.parameters,'rn')
    rn = 1;
else
    rn = settings.parameters.rn;
end
inv_rn = rn\eye(size(rn,1));

%%
% The original publication allows to downsample the data to reduce
% processing time.
%
% Factor is set to 1 by default -> No downsampling.
%
if ~isfield(settings.parameters,'st')
    st = min(1,n);
else
    st = settings.parameters.st;
end
nsmall = ceil(n./st);
wsmall = complex(zeros(settings.parameters.nCh,prod(nsmall),'single'));

%%
% In this section the coil with the highest signal is estimated. Which is
% used in equation [6] of the publication. It is used to calculate the
% field map phase for the j-th coil at the specified location in the FOV.
%
% In case of very low SNR image, we came to the conclusion to smooth the
% image prior to this estimation. This should reduce noise and, therefore,
% reduce the effect of outliers. Here we used a simply 2D Gaussian filter
% with a sigma of 1. Higher sigma results in more blurring of the image.
% There is actually no visual difference in the image, if smoothing is
% applied or not. Therefore, it is currently unclear, if there is any
% improvement (or even a degradation in quality).
%
% Smoothing will be applied to estimate the coil with strongest signal. It
% will not be used for further processing.
%
% If coilComp is enabled maxchannel = 1 will be used. It has to be tested, whtether
% using the coil with the strongest signal will have a beneficial impact.
%
% sigma       = 0.5;
% smoothedIm  = imgaussfilt(real(im),sigma) + 1i * imgaussfilt(imag(im),sigma);

% It is suggested to use the same coil for all slices, because of variation
% of intensity between slices, it helps to get the same phase position. But
% it's not the best solution.
if settings.options.xspace.coilCompression
    maxchannel = 1;
else
    if isfield(settings.parameters,'maxchannel')
        maxchannel = settings.parameters.maxchannel;
    else
        maxchannel = 1;
    end
end

%%
% Adaptive combine uses a sliding window to estimate a weighting factor.
% Referred to block-by-block estimation in the publication. Default in the
% original publication was 4x4.
%
% Bigger sliding window results in better suppression of backgroud noise.
% Too big window size seems to suppress structures as well (64x64 seems to
% be okay for an image matrix of 880x880, 110x110 suppresses structures
% surrounded by mostly noise, e.g. vessels (circle of willis), folded over
% nose tip, etc.
%
%   ->  Signal needs to be strong enough for suppression of background
%       noise. Especially in foot direction.
%
% If no blocksize is specified it will be set to 8x8x3 by default.
%
% If 3D data is input a 3D-sliding will be used.
%
if ~isfield(settings.parameters.xspace,'blockSize')
    bs = min(8,n);
    if n(3) > 1
        bs(3) = min(3,n(3));
    end
else
    bs = settings.parameters.xspace.blockSize;
end

if strcmp(settings.parameters.xspace.combineType, 'bySlice') || n(3) == 1
    bs(3) = 1;
end

%%
% Calulation of the weighting factors.
for z=1:nsmall(3)
    if strcmp(settings.parameters.xspace.combineType, 'bySlice')
        if settings.options.xspace.coilCompression
            tmp = reshape(im(:,:,:,z),settings.parameters.nCh,[]);
            [~,~,V]=svd(tmp*tmp');
            V = V(:,1:nc_svd);
        else
            V = 1;
        end
    elseif ~strcmp(settings.parameters.xspace.combineType, 'bySlice') && z == 1
        if settings.options.xspace.coilCompression
            [~,~,V]=svd(im(:,:)*im(:,:)');
            V = V(:,1:nc_svd);
        else
            V = 1;
        end
    end
    
    for y=1:nsmall(2)
        
        for x=1:nsmall(1)
            cnt = cnt+1;
            ix = [x,y,z];
            imin = max(st.*ix-floor(bs/2)  , 1);
            imax = min(st.*ix+ceil (bs/2)-1, n);
            m1 = V'*reshape(im(:,imin(1):imax(1),imin(2):imax(2),imin(3):imax(3)),settings.parameters.nCh,[]);
            m  = m1*m1';              % Calculate signal covariance
            [v,d]   = eig(inv_rn*m);  % Eigenvector with max eigenval gives
            [~,ind] = max(diag(d));   % the correct combination coeffs.
            tmp = v(:,ind);
            tmp = V*tmp; % transform back to original coil space
            % Correct phase based on coil with max intensity
            tmp = tmp*exp(-1j*angle(tmp(maxchannel)));
            wsmall(:,cnt) = conj(tmp)/(tmp'*inv_rn*tmp);
        end
    end
end
wsmall = reshape(wsmall,[settings.parameters.nCh,nsmall]);

%%
if (settings.parameters.donorm || nargout > 2)
    % This is the normalization proposed in the abstract
    normsmall = reshape(sum(abs(wsmall)).^2,nsmall);
end

%%
% If downsampling enabled above, here the image will be resampled to its
% original matrix.
if (prod(st) ~= 1)
    % Now have to interpolate these weights up to the full resolution. This
    % is done separately for magnitude and phase in order to avoid 0
    % magnitude pixels between +1 and -1 pixels.
    [x,y,z] = ndgrid(1:(nsmall(1)-1)/(n(1)-1):nsmall(1),1:(nsmall(2)-1)/(n(2)-1):nsmall(2),1:(nsmall(3)-1)/(n(3)-1):nsmall(3));
    wfull   = complex(zeros(settings.parameters.nCh, n(1), n(2), n(3),'single'));
    % permute wsmall for faster access
    wsmall  = permute(wsmall,[2 3 4 1]);
    for c=1:settings.parameters.nCh
        if n(3) > 1
            wfull(c,:,:,:) = interpn(abs(wsmall(:,:,:,c)),x,y,z,'linear') .* exp(1j.*interpn(angle(wsmall(:,:,:,c)),x,y,z,'nearest'));
        else % for 2D scans
            wfull(c,:,:,:) = interpn(abs(wsmall(:,:,:,c)),x,y,'linear') .* exp(1j.*interpn(angle(wsmall(:,:,:,c)),x,y,'nearest'));
        end
    end
    if (settings.parameters.donorm || nargout > 2)
        if nz > 1
            norm = interpn(normsmall,x,y,z,'linear');
        else % for 2D scans
            norm = interpn(normsmall,x,y,'linear');
        end
    end
    clear x y z;
else
    wfull  = wsmall;
    if (settings.parameters.donorm || nargout > 2)
        norm  = normsmall;
    end
end
clear wsmall normsmall;

%%
% Applying the weighting factors to the image and combining the channels.
recon = squeeze(sum(wfull.*im));

%%
% If enabled above, here the image is normalized, based on the approach
% given in the publication.
if settings.parameters.donorm
    recon = recon.*norm;
end

% You should carefully read the following terms and conditions before installing or using the
% software. Unless you have entered into a separate written license agreement with
% Universit�t W�rzburg providing otherwise, installation or use of the software indicates your
% agreement to be bound by these terms and conditions.
%
% Use of the software provided with this agreement constitutes your acceptance of these terms.
% If you do NOT agree to the terms of this agreement, promptly remove the software together
% with all copies from your computer. User's use of this software is conditioned upon compliance
% by user with the terms of this agreement.
%
% Upon ordering, downloading, copying, installing or unencrypting any version of the software, you
% are reaffirming that you agree to be bound by the terms of this agreement.
%
% License to use
%
% Universit�t W�rzburg grants to you a limited, non-exclusive, non-transferable and non-assignable
% license to install and use this software for research purposes. Use of this software for any
% diagnostic imaging procedure is strictly forbidden.
%
% License to distribute
%
% Please feel free to offer the non-commercial version of this software on any website, CD, or
% bulletin board, demonstrate the non-commercial version of the software and its capabilities, or
% give copies of the non-commercial version of the software to other potential users, so that others
% may have the opportunity to obtain a copy for use in accordance with the license terms contained
% here.
%
% You agree you will only copy the non-commercial version of the software in whole with this
% license and all delivered files, but not in part.
%
% Termination
%
% This license is effective until terminated. You may terminate it at any point by destroying
% the software together with all copies of the software.
%
% If you have acquired a non-commercial version, the license granted herein shall automatically
% terminate if you fail to comply with any term or condition of this Agreement.
%
% Also, Universit�t W�rzburg has the option to terminate any license granted herein if you fail
% to comply with any term or condition of this Agreement.
%
% You agree upon such termination to destroy the software together with all copies of the software.
%
%
% Copyright
%
% The software is protected by copyright law. You acknowledge that no title to the intellectual
% property in the software is transferred to you. You further acknowledge that title and full
% ownership rights to the software will remain the exclusive property of Universit�t W�rzburg,
% and you will not acquire any rights to the software except as expressly set forth in this
% license. You agree that any copies of the software will contain the same proprietary notices
% which appear on and in the software.
%
% Rent, lease, loan
%
% You may NOT rent, lease or loan the software without first negotiating a specific license
% for that purpose with Universit�t W�rzburg.
%
% No warranties
%
% Universit�t W�rzburg does NOT warrant that the software is error free. Universit�t W�rzburg
% disclaims all warranties with respect to the software, either express or implied, including
% but not limited to implied warranties of merchantability, fitness for a particular purpose and
% noninfringement of third party rights. The software is provided "AS IS."
%
% No liability for consequential damages
%
% In no event will Universit�t W�rzburg be liable for any loss of profits, business, use, or data
% or for any consequential, special, incidental or indirect damages of any kind arising out of
% the delivery or performance or as a result of using or modifying the software, even if
% Universit�t W�rzburg has been advised of the possibility of such damages. In no event will
% Universit�t W�rzburg's liability for any claim, whether in contract, negligence, tort or any
% other theory of liability, exceed the license fee paid by you, if any.
% The licensed software is not designed for use in high-risk activities requiring fail-safe
% performance. Universit�t W�rzburg disclaims any express or implied warranty of fitness for
% high-risk activities.
%
% Severability
%
% In the event of invalidity of any provision of this license, the parties agree that such
% invalidity shall not affect the validity of the remaining portions of this license.
%
