function kspace=applyDenoising(settings, kspace)
% Function to apply various different denoising algorithms.
%
% ************************************************************************
% Version 1.0                                                  15.02.2019
% Falk Luesebrink              falk dot luesebrink at med dot ovgu dot de
% ************************************************************************

if settings.options.kspace.denoise
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Denoise data in image domain before applying any other form of filtering
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Inverse fourier transform k space data into image space
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    data=ifftshift(ifftshift(ifftshift(ifft(ifft(ifft(fftshift(fftshift(fftshift(kspace,1),2),3),[],1),[],2),[],3),1),2),3);
    clear kspace
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Choose denoising algorithm to be used
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if strcmp(settings.options.kspace.filter, 'Wavelet')
        data = ap_denoise_3d(settings, data);
    elseif strcmp(settings.options.kspace.filter, 'BM4D')
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Split complex data into real and imaginary part to denoise separately
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        data_real = real(data);
        data_imag = imag(data);
        clear data
    
        data_real = bm4d(data_real, settings);
        data_imag = bm4d(data_imag, settings);
        data=complex(data_real,data_imag);
%         data = bm4d(data, settings.parameters.denoise.distribution, settings.parameters.denoise.sigma, settings.parameters.denoise.profile, settings.parameters.denoise.do_wiener, settings.parameters.denoise.verbose, settings.parameters.denoise.threshold);
    elseif strcmp(settings.options.kspace.filter, 'Net')
        if verLessThan('matlab', '9.2')
            fprintf('| ***** WARNING ***** You need to have MATLAB 2017b or above to have access to the neural network toolbox.\n');
            fprintf('| ***** WARNING ***** Skipping denoising by neural net.\n');
            return;
        else
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Split complex data into real and imaginary part to denoise separately
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            data_real = real(data);
            data_imag = imag(data);
            clear data
        
            net = denoisingNetwork('DnCNN');

            % The neural net cannot handle data with too low intensity, therefore, it is scaled here.
            % Any better solution? Normalizing each component is not suitable, I suppose.
            % Is there is mathematical correct way to normalize complex valued data?
			% Use scalefac, which is also used for phase unwrapping
            max_real = max(data_real(:));
            max_imag = max(data_imag(:));

            data_real = data_real ./ max_real;
            data_imag = data_imag ./ max_imag;

            for Slice = 1:settings.parameters.pSlc
                data_real(:,:,Slice) = denoiseImage(data_real(:,:,Slice), net);
                data_imag(:,:,Slice) = denoiseImage(data_imag(:,:,Slice), net);
            end

            data_real = data_real .* max_real;
            data_imag = data_imag .* max_imag;

            data=complex(data_real,data_imag);
        end
    elseif strcmp(settings.options.kspace.filter, 'Coupe')
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Split complex data into real and imaginary part to denoise separately
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        data_real = real(data);
        data_imag = imag(data);
        clear data
        
        % Parameter description taken from:
        % sites.google.com/site/pierrickcoupe/CANDLE-User-Manual.pdf
        
        % Smoothing parameter: This parameter controls the amount of denoising
        % that you want to apply. In practice, values between 0.1 and 0.4 are
        % fine for visualization while higher values might be useful for
        % segmentation or registration purposes.
        %     beta = 0.1;
        
        % For large view of fine structures a radius of 1 voxel (ie. patch of
        % 3x3x3 voxels) is recommended while for limited field of view of large
        % structures a radius of 2 voxels is preferred (ie. patch of 5x5x5).
        %     patchradius = 1;
        
        % The radius of the search volume is also a user-controllable option.
        % However, we do not recommend increasing it, as it will produce much
        % longer computational time with minimal denoising improvement.
        %     searchradius = 3;
        
        % Assume Rician (=1) or Gaussian (=0) noise distribution.
        if strcmp(settings.parameters.denoise.distribution, 'Gauss')
            rician = 0;
        else
            rician = 1;
        end
        
        method = settings.parameters.denoise.type;
        
        % fixed range
        map = isnan(data_real(:));
        data_real(map) = 0;
        map = isinf(data_real(:));
        data_real(map) = 0;
        mini = min(data_real(:));
        data_real = data_real - mini;
        maxi = max(data_real(:));
        data_real = data_real*256/maxi;
        
        map = isnan(data_imag(:));
        data_imag(map) = 0;
        map = isinf(data_imag(:));
        data_imag(map) = 0;
        mini = min(data_imag(:));
        data_imag = data_imag - mini;
        maxi = max(data_imag(:));
        data_imag = data_imag*256/maxi;
        
        % Noise estimation
        if method ~= 1
            [hfinal, ho, SNRo, hbg, SNRbg] = MRINoiseEstimation(data_real, rician, settings.parameters.denoise.verbose);
            
            if(isnan(hfinal))
                disp('error during noise estimation')
            end
        end
        
        % Denoising
        if(method==1)
            data_real = MRIDenoisingAONLM(data_real, settings.parameters.denoise.patchradius,  settings.parameters.denoise.searchradius, settings.parameters.denoise.beta, rician, settings.parameters.denoise.verbose);
        end
        
        if(method==2)
            data_real = MRIDenoisingMRONLM(data_real, hfinal, settings.parameters.denoise.beta, settings.parameters.denoise.patchradius,  settings.parameters.denoise.searchradius, rician, settings.parameters.denoise.verbose);
        end
        
        if(method==3)
            data_real = MRIDenoisingONLM(data_real, hfinal, settings.parameters.denoise.beta, settings.parameters.denoise.patchradius, settings.parameters.denoise.searchradius, rician, settings.parameters.denoise.verbose);
        end
        
        if(method==4)
            data_real = MRIDenoisingODCT(data_real, hfinal, settings.parameters.denoise.beta, rician, settings.parameters.denoise.verbose);
        end
        
        if(method==5)
            data_real = MRIDenoisingPRINLM(data_real, hfinal, settings.parameters.denoise.beta, rician, settings.parameters.denoise.verbose);
        end
        map = find(data_real<0);
        data_real(map) = 0;
        
        % Original intensity range
        data_real = data_real * maxi/256;
        data_real = data_real + mini;
        
        % Noise estimation
        if method ~= 1
            [hfinal, ho, SNRo, hbg, SNRbg] = MRINoiseEstimation(data_imag, rician, settings.parameters.denoise.verbose);
            
            if(isnan(hfinal))
                disp('error during noise estimation')
            end
        end
        
        % Denoising
        if(method==1)
            data_imag = MRIDenoisingAONLM(data_imag, settings.parameters.denoise.patchradius,  settings.parameters.denoise.searchradius, settings.parameters.denoise.beta, rician, settings.parameters.denoise.verbose);
        end
        
        if(method==2)
            data_imag = MRIDenoisingMRONLM(data_imag, hfinal, settings.parameters.denoise.beta, settings.parameters.denoise.patchradius,  settings.parameters.denoise.searchradius, rician, settings.parameters.denoise.verbose);
        end
        
        if(method==3)
            data_imag = MRIDenoisingONLM(data_imag, hfinal, settings.parameters.denoise.beta, settings.parameters.denoise.patchradius, settings.parameters.denoise.searchradius, rician, settings.parameters.denoise.verbose);
        end
        
        if(method==4)
            data_imag = MRIDenoisingODCT(data_imag, hfinal, settings.parameters.denoise.beta, rician, settings.parameters.denoise.verbose);
        end
        
        if(method==5)
            data_imag = MRIDenoisingPRINLM(data_imag, hfinal, settings.parameters.denoise.beta, rician, settings.parameters.denoise.verbose);
        end
        map = find(data_imag<0);
        data_imag(map) = 0;
        
        % Original intensity range
        data_imag = data_imag * maxi/256;
        data_imag = data_imag + mini;
        
        data=complex(data_real,data_imag);
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Fourier transform denoised image data back into k space
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    kspace=fftshift(fftshift(fftshift(fft(fft(fft(ifftshift(ifftshift(ifftshift(data,1),2),3),[],1),[],2),[],3),1),2),3);
end