clear; clc; close all;
% Reconstruction pipeline for 3D GRE Siemens raw data acquired with multi
% channel coils.
%
% ************************************************************************
% DISCLAIMER
% 
% PLEASE RESPECT THE COPYWRITE OF INCLUDED THIRD-PARTY TOOLS! THANKS!
% ************************************************************************
% 
% ************************************************************************
% INSTRUCTIONS
%
% Run 'installReconstructionPipeline.m' first and choose a folder where
% data shall be saved to. Then run this script to reconstruct Siemens raw
% data.
% 
% Changes to the default processing can be done either by loading a
% different preset or by changing the settings manually. Presets are stored
% in ./presets/ and to manually change settings open 'setOptions.m' and
% 'setParameters.m', respectively. New presets can be generated using
% 'createPreset.m' using the currently specified options and parameters in
% 'setOptions.m' and 'setParameters.m'.
%
% To process a noise dataset (e.g. for decorrelation of channels), use the
% preset 'noise'. This will change the naming scheme accordingly and run
% the processing pipeline for as long as needed only.
% 
% If you run into errors or anything is unclear, feel free to get in
% contact with me using the address below. Any kind of feedback is highly
% appreciated.
%
% ************************************************************************
% FUNCTIONALITIES
% [*] Reconstruction of raw data files that are larger than memory by
%     writing/reading intermediate results to/from disk.
% [*] Reconstruction of data in parallel by splitting data into chunks.
%     Depending on the processing stage data is split across channels or
%     slices. Number of CPUs and chunks to be processed can be user
%     specified or detected automatically.
% [*] Channel decorrelation by estimating noise covariance matrix (requires
%     noise dataset from the same session).
% [*] Channel combination by sum of squares and adaptive combination.
% [*] Channel combination within adaptive combine slice by slice or in 3D.
% [*] Channel compression within adaptive combine.
% [*] Writing of magnitude, phase, and complex valued NIfTI files.
% [*] 2D GRAPPA reconstruction.
% [*] Application of 3D Tukey filter in k space to reduce Gibbs ringing
%     artifact.
% [*] Application of 3D distortion correction (requires gradient coeffient
%     file of Siemens of your specific gradient system).
% [*] Application of 2D phase unwrapping (per slice).
% [*] Application of denoising in complex domain during reconstruction
%     (includes denoising by BM4D, all NLM denoiser's of Coup� et al. as
%     well as using MATLAB's implementation of DnCNNs utilizing the neural
%     network toolbox).
% [*] Application of bias field correction, segmentation, and creation of a
%     brainmask.
%
% BUGS
% [*] 3D phase unwrapping and wavelet denoising are not working properly.
% [*] If distortion correction is used with gzipped NIfTIs under windows an
%     error occurs (probably use different tool to load data?).
% [*] installReconstructionPipeline has to be run from root directory of
%     the pipeline.
%
% OUTLOOK
% [*] Add version of processing pipeline to settings (and compare version
%     in getStarted or so). 
% [*] Have a look at 'initializeParallel.m' and clean it up. Get nChunks
%     and stepSize, calculate slice and channel memory even no autosetup
%     prepare for output.        
% [*] Rework processing status bar in case of -nodisplay. 
% [*] Change maxMemory parameters relative maximum available memory instead
%     of sabsolute.
% [*] Add safety margin for memory estimations (e.g. 10% of estimation) and
%     add flag to ignore the safefy margin.
% [*] Sanity checks if ~autoSetup for setting up workers, stepSize and
%     nChunks.
% [*] Save output in command window along settings file using 'diary'.
% [*] Implement removal of slice (and phase) oversampling. 
% [*] Display what functions are run in preproc, proc and combine.
% [*] Processing of noise data should be integrated into the processing
%     pipeline without the need of a special preset.
% [*] Write data with correct origin and orientation (currently requires a
%     reference dataset for the correct origin) -> important for distortion
%     correction.
% [*] Add sanity check whether data fits memory for procMemory
% [*] Make use of 'onCleanUp' to close parallel pool if ctrl+c pressed.
% [*] Output BIDS formated.
% [*] If an option is not included in a preset, make sure it is used with a
%     default value (mostly necessary for older versions of a preset?)
%     -> Sanity check during 'setupPipeline.m'
% [*] Make use of BART
% [*] Read in ISMRM rawdata format (and H5?).
% [*] Implement NUFFT for radial trajectories.
% [*] Implement SENSE reconstruction.
% [*] Implement reconstruction in units of SNR.
% [*] Synthesize missing k space data by POCS or homodyne. Code is
%     available, just needs to be reintegrated into the pipeline.
% [*] Auto-adjust filter strength of Tukey filter depending on resolution.
% [*] Create a simple GUI.
% [*] Implement QA (e.g., mriqc.readthedocs.io/en/latest/iqsm/t1w.html)
% [*] Use niftiread and niftiwrite to read and write NIfTI files
%     (introduced in MATLAB 2017b).
% [*] Move from if then else formalism to switches to make the code
%     cleaner.
%
% ************************************************************************
% Version 0.991                                                10.07.2019
% Falk Luesebrink              falk dot luesebrink at med dot ovgu dot de
% ************************************************************************

% Get Version number
versionNumber = '1';

% Use preset settings
preset.usePreset = true;
preset.file      = 'default';

% Process data
name = getStarted;

for numFiles = 1:length(name.FilenameS)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Setup
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    [name, settings, twix_obj]=setupPipeline(name, numFiles, preset, versionNumber);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Retrieve k space from raw data and process into x space
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    [data,settings]=preprocRawdata(name, settings, twix_obj);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Check if noise data has been processed.
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    isNoise = checkNoise(preset);
    if isNoise
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Clean up
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        cleanUp(name, settings, numFiles)
        return;
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Process k space into x space
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    [data,settings]=procRawdata(name, settings, data);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Check if data shall not be combined.
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if settings.options.noCombine
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Clean up
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        cleanUp(name, settings, numFiles)
        return;
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Combine data from x space into image space
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    [data,settings]=procXspace(name, settings, data);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Write NIfTIs
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    writeNIfTI(name, settings, data);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Apply distortion correction (not implemented currently)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    applyDistortionCorrection(name, settings);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Bias correction, segmentation, and brain mask creation
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    applyBiasCorrection(name, settings);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Clean up
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    cleanUp(name, settings, numFiles)
end
clear
