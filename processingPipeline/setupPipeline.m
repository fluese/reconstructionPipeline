function [name, settings, twix_obj]=setupPipeline(name, numFiles, preset, versionNumber)
% Setting up processing pipeline. Using preset or manually specified
% settings.
%
% ************************************************************************
% Version 1.0                                                  15.02.2019
% Falk Luesebrink              falk dot luesebrink at med dot ovgu dot de
% ************************************************************************

fprintf('===============================================================\n');
fprintf('| Reconstruction process started at %s\n', datestr(datetime('now')));
fprintf('===============================================================\n');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Setup options and parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if preset.usePreset
    % Set settings according to m-files
    settings = setOptions;
    settings = setParameters(settings);
    
    % Load settings from preset file
    fpath = fileparts(mfilename('fullpath'));
    tmp = load([fpath '/presets/' preset.file '.mat']);

    % Compare parameters between both settings and add missing parameters
    % in preset to still run pipeline properly.
    f1=fieldnames(settings.parameters);
    f2=fieldnames(tmp.settings.parameters);
    f=intersect(f1,f2);
    tmp=rmfield(settings.parameters, f);
    params=fieldnames(tmp);
    for i=1:length(params)
	    value = getfield(settings.parameters, params{i});
	    tmp.settings.parameters.(params{i}) = value;
    end
    settings.parameters = tmp.settings.parameters;

    % Compare options between both settings and add missing options 
    % in preset to still run pipeline properly. 
    f1=fieldnames(settings.options);
    f2=fieldnames(tmp.settings.options);
    f=intersect(f1,f2);
    tmp=rmfield(settings.options, f);
    opts=fieldnames(tmp);
    for i=1:length(opts)
	    value = getfield(settings.options, opts{i});
	    tmp.settings.options.(opts{i}) = value;
    end
    settings.options = tmp.settings.options;

    % Here general information should be compared. Is it a preset which
    % was created from a previous reconstruction? Do the information
    % match? Should they?

    % Compare version numbers.
else
    settings = setOptions;
    settings = setParameters(settings);
end

settings.versionNumber = versionNumber;
settings.startTime = datestr(datetime('now'));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Setup path for writing files to and path of tissue probability map
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
settings.parameters.path      = '/pool/falk/reconstruction/';                  % Set path were data shall be written in.
settings.parameters.post.path = '/home/luesebri/software/spm12/tpm/TPM.nii';       % Set path to tissue probability model

settings.parameters.numFiles  = numFiles;
name.Filename                 = name.FilenameS{settings.parameters.numFiles,1};
fprintf('| Siemens raw data file to be reconstructed: %s\n', name.Filename)
fprintf('| \n');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Retrieve information of raw data header
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf('| Initializing reconstruction pipeline by parsing header information...\n');
twix_obj = mapVBVD([name.Filepath name.Filename],'removeOS');
fprintf('===============================================================\n\n');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Parameters for reference NIfTI [Currently needed for correct orientation of output file]
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if settings.options.output.refData
    %[refName, refPath] = uigetfile('*.nii;*.nii.gz', 'Select reference data.','MultiSelect', 'off');
    %settings.parameters.output.ref_nifti = [refPath '/' refName];
    
    refs      = cell(8,1);
    refs{1,1} = 'yv98_3512-168_Scanner.nii.gz';
    refs{2,1} = 'yv98_3555-74_Scanner.nii.gz';
    refs{3,1} = 'yv98_3555-75_Scanner.nii.gz';
    refs{4,1} = 'yv98_3589-40_Scanner.nii.gz';
    refs{5,1} = 'yv98_3637-24_Scanner.nii.gz';
    refs{6,1} = 'yv98_3637-25_Scanner.nii.gz';
    refs{7,1} = 'yv98_3681-31_Scanner.nii.gz';
    refs{8,1} = 'yv98_3681-32_Scanner.nii.gz';
    %refs{1,1} = 'sub-01_ses-01_run-01_T1w_bias.nii.gz';
    %refs{2,1} = 'sub-01_ses-02_run-01_T1w_bias.nii.gz';
    %refs{3,1} = 'sub-01_ses-02_run-02_T1w_bias.nii.gz';
    %refs{4,1} = 'sub-01_ses-03_run-01_T1w_bias.nii.gz';
    %refs{5,1} = 'sub-01_ses-04_run-01_T1w_bias.nii.gz';
    %refs{6,1} = 'sub-01_ses-04_run-02_T1w_bias.nii.gz';
    %refs{7,1} = 'sub-01_ses-05_run-01_T1w_bias.nii.gz';
    %refs{8,1} = 'sub-01_ses-05_run-02_T1w_bias.nii.gz';

    refName = refs{numFiles,1};
    refPath = '/pool/public/data/250um/';
    %refPath = '/pool/falk/data/RegistrationFinal/data/';

    if refName == 0
        fprintf('| No reference file selected. Contiuing without reference.\n')
        settings.options.output.refData = false;
    else
        [~,~,ext]=fileparts(refName);
        if strcmp(ext, '.gz')   
            gunzip([refPath '/' refName]);
            refName=refName(1:end-3);
        end
        settings.parameters.output.ref_nifti = [refPath '/' refName];
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Setup parameters and name based on MHD
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
settings.parameters.nCh  = twix_obj.image.NCha;
settings.parameters.nRO  = twix_obj.hdr.Meas.NImageCols;
settings.parameters.nPE  = twix_obj.hdr.Meas.NImageLins;
settings.parameters.nSlc = twix_obj.hdr.Meas.Partitions;
name.nameVol             = twix_obj.hdr.Meas.tPatientName;
name.nameSession         = ['-' num2str(twix_obj.hdr.Config.MeasUID)];

% Parameters used in processing partial Fourier data
settings.parameters.pRO  = floor(twix_obj.image.NCol/2);                        % Length of partially sampled points in readout direction (without RO oversampling)
settings.parameters.pPE  = twix_obj.image.NLin;                                 % Length of partially sampled points in phase encoding direction
settings.parameters.pSlc = twix_obj.image.NPar;                                 % Length of partially sampled points in slice direction
settings.parameters.uRO  = settings.parameters.nRO  - settings.parameters.pRO;  % Length of undersampled points in readout direction
settings.parameters.uPE  = settings.parameters.nPE  - settings.parameters.pPE;  % Length of undersampled points in phase encoding direction
settings.parameters.uSlc = settings.parameters.nSlc - settings.parameters.pSlc; % Length of undersampled points in slice direction

% Setup partial Fourier
if settings.parameters.uRO ~= 0 || settings.parameters.uPE ~= 0 || settings.parameters.uSlc ~= 0
    settings.options.kspace.partialFourier    = true;
    settings.parameters.kspace.partialFourier = 'zero';
elseif settings.parameters.uRO == 0 && settings.parameters.uPE == 0 && settings.parameters.uSlc == 0
    settings.options.kspace.partialFourier    = false;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Create directory where data is stored in
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
name.namePath = [settings.parameters.path name.nameVol];
if ~exist(name.namePath, 'dir')
    mkdir(name.namePath)
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Get orientation information
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
settings = getOrientation(twix_obj, settings);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Is GRAPPA?
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if isfield(twix_obj, 'refscan')
    settings.options.kspace.GRAPPA = true;
else
    settings.options.kspace.GRAPPA = false;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Channel decorrelation (optional)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% if isfield(twix_obj, 'noise')
%     settings.options.kspace.decorrelateChannels    = true;
%     settings.parameters.kspace.decorrelateChannels = 'noiseData';
% end
settings = decorrelationMatrix(name, settings, twix_obj);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Sanity checks (Should be put together in a function once more checks are
% necessary. Also each check should be commented!)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% if settings.options.output.writeChannelwise && settings.options.procMemory
%     settings.options.saveIntermediate = true;
%     settings.options.force            = true;
% end

% In case proc memory is set to false, so should be the saveIntermediate option.
if ~settings.options.procMemory
	settings.options.saveIntermediate = false;
end

if settings.options.saveIntermediate && settings.options.procMemory
    settings.options.keepFiles = true;
elseif settings.options.procMemory && ~settings.options.saveIntermediate
    settings.options.keepFiles = false;
end

% Change the number of CPUs used during each stage.
if settings.parameters.maxCPUs ~= Inf && ~settings.options.autoSetup
    settings.parameters.preprocWorkers = settings.parameters.maxCPUs;
    settings.parameters.combineWorkers = settings.parameters.maxCPUs;
    settings.parameters.procWorkers    = settings.parameters.maxCPUs;
    
    if settings.parameters.maxCPUs == 1
        settings.options.procParallel  = false;
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Setup naming scheme
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
name = setName(name, settings, preset);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Display options
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
displaySettings(name, settings, preset);

