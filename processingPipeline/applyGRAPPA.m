function res=applyGRAPPA(kData,kCalib,kSize,lambda)
% Function to run GRAPPA reconstruction in reconstruction pipeline.
% Original code copyrighted by Michael Lustig. Edited by Hendrik Mattern
% and Falk Luesebrink.
%
% **************************************************
% 
% res = GRAPPA(kData,kCalib,kSize,lambda [, disp)
%
% This is a GRAPPA reconstruction algorithm that supports 
% arbitrary Cartesian sampling. However, the implementation
% is highly inefficient in Matlab because it uses for loops. 
% This implementation is very similar to the GE ARC implementation.
%
% The reconstruction looks at a neighborhood of a point and
% does a calibration according to the neighborhood to synthesize
% the missing point. This is a k-space varying interpolation.
% A sampling configuration is stored in a list, and retrieved
% when needed to accelerate the reconstruction (a bit)
%
% Inputs: 
%       kData     - [Size x, Size y, num coils] 2D multi-coil k-space data to reconstruct from.
%                   Make sure that the missing entries have exact zeros in them.
%       kCalib    - calibration data (fully sampled k-space)
%       kSize     - size of the 2D GRAPPA kernel [kx, ky]
%       lambda    - Tykhonov regularization for the kernel calibration.
%       dispp      - Figure number to display images as they are
%                   reconstructed
% Outputs:
%       res       - k-space data where missing entries have been filled in.
%
% Example:
%   [x,y] = meshgrid(linspace(0,1,128));
%   % Generate fake Sensitivity maps
%   sMaps = cat(3,x.^2,1-x.^2,y.^2,1-y.^2);
%   % generate 4 coil phantom
%   imgs = repmat(phantom(128),[1,1,4]).*sMaps;
%   DATA = fft2c(imgs);
%   % crop 20x20 window from the center of k-space for calibration
%   kCalib = crop(DATA,[20,20,4]);
%
%   %calibrate a kernel
%   kSize = [5,5];
%   coils = 4;
%   
%   % undersample by a factor of 2
%   DATA(1:2:end,2:2:end,:) = 0;
%   DATA(2:2:end,1:2:end,:) = 0;
%   
%   %reconstruct:
%   [res] = GRAPPA(DATA,kCalib, kSize, 0.01);
%   figure, imshow(cat(2,sos(imgs), 2*sos(ifft2c(DATA)), sos(ifft2c(res))),[]);
%   title('full,  zero-fill,   result')
%   
%       
%
% (c) Michael Lustig 2008
%
% **************************************************
%
% Edited by Falk Luesebrink 20.02.2019 (falk dot luesebrink at med dot ovgu dot de) 


szData  = size(kData); % get sizes
szCalib = size(kCalib);
res = kData*0;

AtA = dat2AtA(kCalib, kSize); % build coil calibrating matrix

par = false;
if par
    parfor n=1:szData(3) 	
        resCoil = ARC(kData, AtA, kSize, n, lambda); % reconstruct single coil image
        res(:,:,n) = resCoil;
    end
else
    for n=1:szData(3) 	
        resCoil = ARC(kData, AtA, kSize, n, lambda); % reconstruct single coil image
        res(:,:,n) = resCoil;
    end    
end

lineIdx=floor((szData(2)-szCalib(2))./2);
res(:, lineIdx+1:lineIdx+szCalib(2),:)=kCalib;


function [res] = ARC(kData, AtA, kSize, c, lambda)
[sx,sy,nCoil] = size(kData);


kData = zpad(kData,[sx+kSize(1)-1, sy+kSize(2)-1, nCoil]);

dummyK = zeros(kSize(1),kSize(2),nCoil); dummyK((end+1)/2,(end+1)/2, c) = 1;
idxy = find(dummyK);

res = zeros(sx,sy);

MaxListLen = 100;
LIST = zeros(kSize(1)*kSize(2)*nCoil,MaxListLen);
KEY =  zeros(kSize(1)*kSize(2)*nCoil,MaxListLen);
count = 0;

for y = 1:sy
	for x=1:sx
		tmp = kData(x:x+kSize(1)-1,y:y+kSize(2)-1,:);
		pat = abs(tmp)>0;
		if pat(idxy) | sum(pat)==0
			res(x,y) = tmp(idxy);
		else
			key = pat(:);
            idx = 0;
			for nn=1:size(KEY,2)
				if sum(key==KEY(:,nn))==length(key)
				   idx = nn;
				   break;
			   	end
			end
			if idx == 0
				count = count + 1;
				kernel = calibrate(AtA,kSize,nCoil,c,lambda,pat);
				KEY(:,mod(count,MaxListLen)+1) = key(:);
				LIST(:,mod(count,MaxListLen)+1) = kernel(:);
			else
				kernel = LIST(:,idx);
			end
			res(x,y) = sum(kernel(:).*tmp(:));
		end

	end
end