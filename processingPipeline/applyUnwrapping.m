function data=applyUnwrapping(settings, data)
% Function to unwrap complex data (code provided by Andre Pampel of 
% Max-Planck-Institute Leipzig)
%
% ************************************************************************
% Version 1.0                                                  15.02.2019
% Falk Luesebrink              falk dot luesebrink at med dot ovgu dot de
% ************************************************************************

if settings.options.kspace.unwrap	
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Inverse fourier transform k space data into image space
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    data=ifftshift(ifftshift(ifftshift(ifft(ifft(ifft(fftshift(fftshift(fftshift(data,1),2),3),[],1),[],2),[],3),1),2),3);
	
    for Slice = 1:settings.parameters.pSlc
        data(:,:,Slice)=tv_pc_silent(data(:,:,Slice),0.03,3,settings.parameters.kspace.scalefac);
    end
    
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% Fourier transform denoised image data back into k space
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    data=fftshift(fftshift(fftshift(fft(fft(fft(ifftshift(ifftshift(ifftshift(data,1),2),3),[],1),[],2),[],3),1),2),3);
end