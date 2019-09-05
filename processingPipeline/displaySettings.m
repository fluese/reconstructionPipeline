function displaySettings(name, settings, preset)
% Function to display settings and image information.
%
% ************************************************************************
% Version 1.0                                                  15.02.2019
% Falk Luesebrink              falk dot luesebrink at med dot ovgu dot de
% ************************************************************************

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Display image information
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf('===============================================================\n');
fprintf('| Image information\n');
fprintf('===============================================================\n');
fprintf('|\n')
fprintf('| Dimensions: %i x %i x %i x %i\n', settings.parameters.nCh, settings.parameters.nRO, settings.parameters.nPE, settings.parameters.nSlc);
if settings.parameters.nRO > settings.parameters.pRO || settings.parameters.nPE > settings.parameters.pPE || settings.parameters.nSlc > settings.parameters.pSlc
	fprintf('| Dimensions (partially acquired): %i x %i x %i x %i\n', settings.parameters.nCh, settings.parameters.pRO, settings.parameters.pPE, settings.parameters.pSlc);
end
fprintf('| Resolution: %.2f x %.2f x %.2f\n', settings.parameters.Resolution(1), settings.parameters.Resolution(2), settings.parameters.Resolution(3));
fprintf('|\n')
fprintf('===============================================================\n\n');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Display settings
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf('===============================================================\n');
fprintf('| Settings\n');
fprintf('===============================================================\n');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Loading preset?
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if preset.usePreset
    fprintf(['| Loading preset: ' preset.file '.\n'])
    fprintf('| \n');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Auto setup?
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if settings.options.autoSetup
    fprintf('| Automatically set maximum number of CPUs and/or slices per chunk ')
    if settings.parameters.maxMemory == Inf
        fprintf('according to available memory.\n');	
    else
        fprintf('limiting used memory to %i GB.\n', settings.parameters.maxMemory);
    end
else
    fprintf('| Manually set maximum number of CPUs and/or slices per chunk.\n')
end

% Memory needs to estimated in case of manual setup.
%fprintf('| \n');
%fprintf('| Number of CPUs for preprocessing: %i\n', settings.parameters.preprocWorkers);
%fprintf('| Number of CPUs for processing: %i\n', settings.parameters.procWorkers);
%fprintf('| Number of CPUs for channel combination: %i\n', settings.parameters.combineWorkers);
%fprintf('| \n');
%fprintf('| Available memory: %.2f\n', settings.parameters.availableMem);
%fprintf('| Estimated slice memory: %.2f\n', settings.parameters.SliceMem*settings.parameters.nChunk);
%fprintf('| Estimated channel memory: %.2f\n', settings.parameters.ChannelMem*seetings.parameters.nChunks);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Processing in parallel?
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if settings.options.procParallel
	fprintf('| Data are processed in parallel.\n')
else
	fprintf('| Data are not processed in parallel.\n')
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Combination of channels
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if settings.options.noCombine
	fprintf('| Channels will not be combined.\n');
elseif strcmp(settings.options.xspace.combineChannels, 'adaptComb')
	if strcmp(settings.parameters.xspace.combineType,'byChunk')
	        fprintf('| Channels are combined by adaptive combine with a block size of %ix%ix%i.\n', settings.parameters.xspace.blockSize);
    else
		fprintf('| Channels are combined by adaptive combine per slice with a block size of %ix%i.\n', settings.parameters.xspace.blockSize);
	end
elseif strcmp(settings.options.xspace.combineChannels, 'SoS')
    fprintf('| Channels are combined by sum of squares.\n');
end

fprintf('|\n');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Decorrelation of channels
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if settings.options.kspace.decorrelateChannels && ~isscalar(settings.parameters.W)
    fprintf('| Estimate decorrelation matrix ');

    if strcmp(settings.parameters.kspace.decorrelateChannels, 'slice')
        fprintf('based on first slice.\n');
    elseif strcmp(settings.parameters.kspace.decorrelateChannels, 'noiseScan')
        fprintf('based on separately acquired noise data.\n');
    elseif strcmp(settings.parameters.kspace.decorrelateChannels, 'noiseData')
        fprintf('based on noise data within scan.\n');
    end
elseif settings.options.kspace.decorrelateChannels && isscalar(settings.parameters.W)
    fprintf('| Option to decorrelated channels was set to true, however, no noise data was found...\n')
    fprintf('| Noise covariance to decorrelate channels will not be estimated.\n');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% GRAPPA reconstruction?
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if settings.options.kspace.GRAPPA
	fprintf('| Running GRAPPA reconstruction.\n');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Denoising?
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if settings.options.kspace.denoise
    fprintf('| Data are denoised using ');

    if strcmp(settings.options.kspace.filter, 'Wavelet')
        fprintf('wavelets.\n')
        
    elseif strcmp(settings.options.kspace.filter, 'BM4D')
		if strcmp(settings.parameters.denoise.profile, 'np')
			fprintf('BM4D with normal profile.\n')
		elseif strcmp(settings.parameters.denoise.profile, 'lc')
			fprintf('BM4D with low complexity profile.\n')
		elseif strcmp(settings.parameters.denoise.profile, 'mp')
			fprintf('BM4D with modified profile.\n')
		elseif strcmp(settings.parameters.denoise.profile, 'manual')
			fprintf('BM4D with manual settings.\n')
		end
			
    elseif strcmp(settings.options.kspace.filter, 'Net')
        fprintf('DnCNN.\n')
        
    elseif strcmp(settings.options.kspace.filter, 'Coupe')
	if settings.parameters.denoise.type == 1
           fprintf('AONLM ')
	elseif settings.parameters.denoise.type == 2
           fprintf('MR-ONLM ')
	elseif settings.parameters.denoise.type == 3
           fprintf('ONLM ')
	elseif settings.parameters.denoise.type == 4
           fprintf('ODCT ')
	elseif settings.parameters.denoise.type == 5
           fprintf('PRINLM ')
	end
	fprintf('with a beta of %.1f.\n', settings.parameters.denoise.beta)
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% K space filtering?
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if settings.options.kspace.tukeyWindow
    fprintf('| Application of Tukey window with a taper length of %.2f.\n', settings.parameters.kspace.tukey);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Partial Fourier reconstruction?
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if settings.options.kspace.partialFourier
    fprintf('| Partial Fourier reconstruction by ')
    if strcmp(settings.parameters.kspace.partialFourier, 'zero')
        fprintf('zero filling.\n')
    elseif strcmp(settings.parameters.kspace.partialFourier, 'conj')
        fprintf('***** Currently not implemented. *****\n')
    elseif strcmp(settings.parameters.kspace.partialFourier, 'pocs')
        fprintf('***** Currently not implemented. *****\n')
    elseif strcmp(settings.parameters.kspace.partialFourier, 'homo')
        fprintf('***** Currently not implemented. *****\n')
    end
end

fprintf('| \n');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Bias field correction (and segmentation)?
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if settings.options.post.biasfield
    if ~settings.parameters.SPM
        fprintf('| SPM is not installed or TPM path is not setup properly. Will not run inhomogeneity correction.\n');
    else
        if settings.options.post.segmentation
            fprintf('| Bias field correction and segmentation will be conducted.\n')
        else
            fprintf('| Bias field correction will be conducted.\n')
        end
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Distortion correction?
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if settings.options.post.distortionCorrection
    fprintf('| Distortion correction will be conducted.\n')
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Files will be written like this and that.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf('|\n');
disp(['| Files will be written in subject folder: ' name.namePath])

if settings.options.output.writeChannelwise
    fprintf('| Unprocessed data will be written to disk as NIfTI files per channel.\n')
end

if settings.options.output.writeNIfTI
    % Could be done better with by settings up an extension string which
    % varies depending on settings.options.output.gzip. Would clear up
    % everything.
    if settings.options.output.writePhase
		if settings.options.output.gzip
			fprintf(['| Magnitude and phase data will be written as: ' name.nameFile '.nii.gz\n'])
		else
			fprintf(['| Magnitude and phase data will be written as: ' name.nameFile '.nii\n'])
		end
	else
		if settings.options.output.gzip
			fprintf(['| Magnitude data will be written as: ' name.nameFile '.nii.gz\n'])
		else
			fprintf(['| Magnitude data will be written as: ' name.nameFile '.nii\n'])
		end
    end
end

if settings.options.output.writeComplex
	if settings.options.output.gzip
		fprintf(['| Complex data will be written as: ' name.nameFile '_complex.nii.gz\n'])
	else
		fprintf(['| Complex data will be written as: ' name.nameFile '_complex.nii\n'])
	end
end

if settings.options.output.refData
	if settings.options.output.gzip
		fprintf(['| Reference data will be used to set the correct origin: ' name.nameFile '_withRef.nii.gz\n'])
	else
		fprintf(['| Reference data will be used to set the correct origin: ' name.nameFile '_withRef.nii\n'])
	end
end

if settings.options.keepFiles
    fprintf('| Intermediately written *.mat files will be kept.\n')
elseif ~settings.options.keepFiles && ~settings.options.procMemory
    fprintf('| Intermediately *.mat files will be deleted during clean up.\n')
end

if settings.options.output.gzip
    fprintf('| All NIfTI files will be gzipped.\n')
end

if settings.options.force 
	fprintf('| ****** WARNING! Forcing to overwrite all existing files.\n')
elseif settings.options.kspace.force
	fprintf('| ****** WARNING! Forcing to overwrite all existing files.\n')
elseif settings.options.xspace.force 
	fprintf('| ****** WARNING! Forcing to overwrite x space mat files and beyond.\n')
elseif settings.options.imspace.force
	fprintf('| ****** WARNING! Forcing to overwrite image space mat files and beyond.\n')
elseif settings.options.NIfTI.force
	if settings.options.output.writeNIfTI && settings.options.output.writePhase
		fprintf('| ****** WARNING! Forcing to overwrite magnitude and phase NIfTI files.\n')
	elseif settings.options.output.writePhase
		fprintf('| ****** WARNING! Forcing to overwrite phase NIfTI files.\n')
	elseif settings.options.output.writeNIfTI
		fprintf('| ****** WARNING! Forcing to overwrite magnitude NIfTI files.\n')
	end
end

fprintf('===============================================================\n\n');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% End
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
