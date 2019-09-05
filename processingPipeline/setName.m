function name=setName(name, settings, preset)
% Function to set up naming scheme.
%
% ************************************************************************
% Version 1.0                                                  15.02.2019
% Falk Luesebrink              falk dot luesebrink at med dot ovgu dot de
% ************************************************************************

% Prefix
name.namePrefix     = '';
name.nameSuffix     = '';

% Denoise data
if settings.options.kspace.denoise
    if strcmp(settings.options.kspace.filter, 'Wavelet')
        name.nameDenoise = '_denoised-Wavelet';
    elseif strcmp(settings.options.kspace.filter, 'BM4D')
        name.nameDenoise = strrep(['_denoised-BM4D_distribution-' settings.parameters.denoise.distribution '_profile-' settings.parameters.denoise.profile], '.', '');
    elseif strcmp(settings.options.kspace.filter, 'Net')
        name.nameDenoise = '_denoised-Net';
    elseif strcmp(settings.options.kspace.filter, 'Coupe')
        if settings.parameters.denoise.type == 1
               type='AONLM';
        elseif settings.parameters.denoise.type == 2
               type='MR-ONLM';
        elseif settings.parameters.denoise.type == 3
               type='ONLM';
        elseif settings.parameters.denoise.type == 4
               type='ODCT';
        elseif settings.parameters.denoise.type == 5
               type='PRINLM';
        end
        name.nameDenoise = strrep(['_denoised-' type '_beta-' num2str(settings.parameters.denoise.beta) '_distribution-' settings.parameters.denoise.distribution], '.', '');
    end
else
    name.nameDenoise = '';
end

% Decorrelate channels
if settings.options.kspace.decorrelateChannels
    if strcmp(settings.parameters.kspace.decorrelateChannels, 'slice')
        name.nameDecor = '_decorrelated-slice';
    elseif strcmp(settings.parameters.kspace.decorrelateChannels, 'noiseData')
        name.nameDecor = '_decorrelated-noiseData';
    elseif strcmp(settings.parameters.kspace.decorrelateChannels, 'noiseScan')
        name.nameDecor = '_decorrelated-noiseScan';
    end
else
    name.nameDecor = '';
end

% Partial Fourier reconstruction
if settings.options.kspace.partialFourier
    if strcmp(settings.parameters.kspace.partialFourier, 'zero')
        name.namePF = ['_PF-zero'];
    elseif strcmp(settings.parameters.kspace.partialFourier, 'conj')
        name.namePF = ['_PF-conj'];
    elseif strcmp(settings.parameters.kspace.partialFourier, 'pocs')
        name.namePF = ['_PF-POCS'];
    elseif strcmp(settings.parameters.kspace.partialFourier, 'homo')
        name.namePF = ['_PF-homo'];
    end
else
    name.namePF = '';
end

% Phase unwrapping
if settings.options.kspace.unwrap
    name.nameUnwrap = '_unwrapped';
else
    name.nameUnwrap = '';
end

% Retrieve channel with maximum intensity for phase correction
if settings.options.xspace.getMaxChannel && strcmp(settings.options.xspace.combineChannels, 'adaptComb')
    name.nameMax = '_maxChannel';
else
    name.nameMax = '';
end

% Application of Tukey filter for removal of Gibbs ringing
if settings.options.kspace.tukeyWindow
    name.nameTukey = strrep(['_Tukey-' num2str(settings.parameters.kspace.tukey)], '.', '');               % Name refering to k space filtering
else
    name.nameTukey   = '';
end

% Type of channel combination
if strcmp(settings.options.xspace.combineChannels, 'adaptComb')
    name.type = '_adaptComb';
    if strcmp(settings.parameters.xspace.combineType, 'bySlice')
        name.nameComb = strrep([name.type '-' num2str(settings.parameters.xspace.blockSize(1)) 'x' num2str(settings.parameters.xspace.blockSize(2))], '.', ''); % Name refering to coil combination
    elseif strcmp(settings.parameters.xspace.combineType, 'byChunk')
        name.nameComb = strrep([name.type '-' num2str(settings.parameters.xspace.blockSize(1)) 'x' num2str(settings.parameters.xspace.blockSize(2)) 'x' num2str(settings.parameters.xspace.blockSize(3))], '.', ''); % Name refering to coil combination
    end
elseif strcmp(settings.options.xspace.combineChannels, 'SoS')
    name.type  = '_SoS';
    name.nameComb = name.type;
end

if settings.options.output.shortName
    name.namePath = [settings.parameters.path name.nameVol];
    name.kspace   = [name.namePrefix 'kspace' name.nameSession name.nameSuffix];
    name.xspace   = [name.namePrefix 'xspace' name.nameSession name.nameSuffix];
    name.imspace  = [name.namePrefix 'imspace' name.nameSession name.nameSuffix];
    name.nameFile = [name.namePrefix name.nameVol name.nameSession name.nameSuffix];
	    if settings.options.output.gzip
	        name.NIfTI    = [name.namePath '/' name.nameFile '.nii.gz'];
	    else
	        name.NIfTI    = [name.namePath '/' name.nameFile '.nii'];
	    end
elseif settings.options.output.shortFilename
	name.namePath = [settings.parameters.path name.nameVol];
        name.kspace   = [name.namePrefix 'kspace' name.nameSession name.nameDecor name.nameSuffix];
        name.xspace   = [name.namePrefix 'xspace' name.nameSession name.nameDecor name.nameDenoise name.nameTukey name.namePF name.nameUnwrap name.nameSuffix];
        name.imspace  = [name.namePrefix 'imspace' name.nameSession name.nameDecor name.nameDenoise name.nameTukey name.namePF name.nameUnwrap name.nameMax name.nameComb name.nameSuffix];
	name.nameFile = [name.namePrefix name.nameVol name.nameSession name.nameSuffix];
else
	name.namePath = [settings.parameters.path name.nameVol];
	name.kspace   = [name.namePrefix 'kspace' name.nameSession name.nameDecor name.nameSuffix];
	name.xspace   = [name.namePrefix 'xspace' name.nameSession name.nameDecor name.nameDenoise name.nameTukey name.namePF name.nameUnwrap name.nameSuffix];
	name.imspace  = [name.namePrefix 'imspace' name.nameSession name.nameDecor name.nameDenoise name.nameTukey name.namePF name.nameUnwrap name.nameMax name.nameComb name.nameSuffix];
	name.nameFile = [name.namePrefix name.nameVol name.nameSession name.nameDecor name.nameDenoise name.nameTukey name.namePF name.nameUnwrap name.nameMax name.nameComb name.nameSuffix];
	if settings.options.output.gzip
	    name.NIfTI    = [name.namePath '/' name.nameFile '.nii.gz'];
	else
	    name.NIfTI    = [name.namePath '/' name.nameFile '.nii'];
	end
end

if preset.usePreset && strcmp(preset.file, 'noise')
    name.kspace     = 'Noise';
    name.xspace     = 'Noise_xspace';
	name.imspace    = 'Noise_imspace';
    name.nameFile   = [name.nameVol '_Noise'];
end
