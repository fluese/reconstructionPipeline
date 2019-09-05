function settings=setOptions
% Function to set up options.
%
% ************************************************************************
% Version 1.0                                                  15.02.2019
% Falk Luesebrink              falk dot luesebrink at med dot ovgu dot de
% ************************************************************************

% Global options
settings.options.procParallel                   = true;              % Reconstruct data and combine channels in parallel. (Mostly included for debugging, as debugging using the parallel toolbox is a pain.)
settings.options.autoSetup                      = true;              % Automatically setup number of workers to process chunks of data in parallel. If false the parameters can be set up in 'setParameters.m'
settings.options.memSafetyMargin                = true;              % Use a safety margin for memory estimation. Default: 90% of available memory. Can be changed in setParameters.m
settings.options.keepFiles                      = false;             % Keep all intermediately created files during reconstruction.
settings.options.procMemory                     = false;             % Instead of reading intermediate steps from disk, keep everything in memory to increase processing speed.
settings.options.saveIntermediate               = false; 			 % If procMemory = true: Write intermediate results to disk.
settings.options.force                          = false;             % Force to process data. All files will be overwritten.
settings.options.noCombine                      = false;             % Stop reconstruction before combination of channels

% Force writing data
settings.options.kspace.force                   = false;             % [Pre-processing kspace] Force overwriting k space data (will set settings.options.xspace.force and settings.options.imspace.force = true)
settings.options.xspace.force                   = false;             % [Processing kspace] Force overwriting x space data (will set settings.options.imspace.force = true)
settings.options.imspace.force                  = false;             % [Channel combination] Force overwriting image space data (will set settings.options.NIfTI.force = true)
settings.options.NIfTI.force                    = false;             % [NIfTI generation] Force overwriting NIfTI file

% How to reconstruct k space data into x space
settings.options.kspace.decorrelateChannels     = false;             % Estimate noise covariance to decorrelate channels (Needs to have processed noise data)
settings.options.kspace.unwrap                  = false;             % Conduct phase unwrapping using total variation (based on code provided by Max Planck Institute in Leipzig)
settings.options.kspace.denoise                 = false;             % Denoise complex data per channel dopen uring reconstruction.
settings.options.kspace.filter                  = 'BM4D';            % Choose your denoising filter (currently implemented 'BM4D', 'Net' and 'Wavelet', 'Coupe' is work in progress.)
settings.options.kspace.tukeyWindow             = true;              % Apply Tukey window on k space data to remove Gibbs ringing.
settings.options.kspace.variableAveraging       = false;             % Sum k-space along averages and divide by number of averages
settings.options.kspace.acquisitionWeighting    = false;             % Multiply k-space by weighting factor (Hann function)
% settings.options.kspace.interpolate             = false;           % Interpolate k-space data (not implemented yet)

% How to combine channels from x space into image space
settings.options.xspace.combineChannels         = 'adaptComb';       % Type of channel combination. ('SoS' = Sum of squares, 'adaptComb' = Adaptive combine.) 
settings.options.xspace.coilCompression         = false;             % Use coil compression in adaptive combine (potentially makes sense to implement much ealier to speed up reconstruction).
settings.options.xspace.getMaxChannel           = true;              % Load entire dataset and retrieve the channel with the highest intensity for phase correction in adaptive combine.

% Writing of output
settings.options.output.writeNIfTI              = true;              % Write magnitude NIfTI (has to be set to true, for all other files to be written)
settings.options.output.writePhase              = false;             % Write phase NIfTI
settings.options.output.writeComplex            = false;             % Write complex valued NIfTI (as a time-series with real part as first time point and imaginary part as second time point)
settings.options.output.gzip                    = true;              % Gzip output files
settings.options.output.refData                 = false;             % Provide reference NIfTI data for correct origin and orientation.
settings.options.output.trueComplex             = false;             % Write true complex valued NIfTI
settings.options.output.writeChannelwise        = false;             % Save unprocessed complex data per channel to disk
settings.options.output.shortName               = true;              % Instead of using parameters in folder and filenames, use shortest possible unambigous name.
settings.options.output.shortFilename           = false;             % Use short name for NIfTI only
% settings.options.output.BIDS                    = false;           % Save data BIDS conform (not implemented yet)

% Post-processing
settings.options.post.distortionCorrection      = false;             % Apply 3D distortion correction
settings.options.post.biasfield                 = true;             % Apply bias field correction of SPM12
settings.options.post.segmentation              = false;             % Write segmentation of GM, WM, and CSF
settings.options.post.mask                      = false;             % Create brainmask (requires parameter.segmentation = true)
settings.options.post.keepSegmentation          = false;             % Keep segmentation or remove it after creation of brainmask
settings.options.post.keepBiascorrection        = false;             % Keep bias corrected file or remove it after creation of brainmask