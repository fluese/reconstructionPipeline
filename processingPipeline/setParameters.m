function settings=setParameters(settings)
% Function to set up parameters. Should be treated like expert options.
%
% ************************************************************************
% Version 1.0                                                  15.02.2019
% Falk Luesebrink              falk dot luesebrink at med dot ovgu dot de
% ************************************************************************

% Global parameters
settings.parameters.maxCPUs                    = Inf;                              % Set the maximum number of CPUs globally.
settings.parameters.preprocWorkers             = Inf;                              % Number of CPUs to process chunks of slices in parallel. In case of settings.options.autoSetup: Inf = Use maximum amount of physical cores.
settings.parameters.combineWorkers             = Inf;                              % Number of CPUs to process chunks of slices in parallel. In case of settings.options.autoSetup: Inf = Use maximum amount of physical cores.
settings.parameters.procWorkers                = Inf;                              % Number of CPUs to process chunks of channels in parallel. In case of settings.options.autoSetup: Inf = Use maximum amount of physical cores.
settings.parameters.maxMemory                  = Inf;                              % Only in case of settings.options.autoSetup = true. Inf = Use memory available otherwise limit the maximum amount of memory used in GB. 
settings.parameters.stepSize                   = [];
settings.parameters.nChunks                    = [];
settings.parameters.safetyMargin               = 1;                                % Limit the maximum amount of available memory given from 0 to 1.
settings.parameters.scalefac                   = 1000;

% Parameter for reconstruction
settings.parameters.kspace.scalefac 	       = 1/1e-11; 							% Scaling factor for phase unwrapping (and potentiatlly denosing)
settings.parameters.kspace.decorrelateChannels = 'noiseScan';                       % Type of decorrelation method. ('noiseScan' = based on separately acquired noise data, 'noiseData' = based on noise data within raw data file)
settings.parameters.kspace.tukey               = 0.05;                               % Taper length of Tukey filter for k space filtering (between 0 and 1 with 0: Rect window, 1: Hanning window)
settings.parameters.kspace.partialFourier      = 'zero';                            % Apply partial Fourier reconstruction. ('zero' = zero filling, 'pocs' = POCS, 'homo' = homodyne, 'conj' = conjugate synthesize.)

% Parameters for channel combination
settings.parameters.xspace.blockSize           = [4 4 4];                           % Block size for adaptive combine
settings.parameters.xspace.combineType         = 'byChunk';                         % Combine channels in in 2D ('bySlice') or 3D ('byChunk') for adaptive combine

% Parameters for bias field correction:
% See 'bias_correction.m' for in-depth explantion of settings.parameters.
settings.parameters.post.fwhm                  = 30;                                % FWHM (default: 60)
settings.parameters.post.reg                   = 0.001;                             % Regularization (default: 0.001)
settings.parameters.post.samp                  = 2;                                 % Sampling distance (default: 3)

% Parameters for brainmask
settings.parameters.post.dilate                = 10;                                % Dilate brainmask by cube with edge length of parameter.dilate
settings.parameters.post.erode                 = 10;                                % If greater than 0, erode brainmask after dilation

% Parameters for 3D phase correction
settings.parameters.pk.name                    = {'gauss','gauss','gauss'};
settings.parameters.pk.width                   = [8 8 8];
settings.parameters.pk.shifts                  = [0 0 0];
settings.parameters.pk.scale                   = false;
settings.parameters.pk.writephasemap           = false;

% General parameters for denoising the data
settings.parameters.denoise.distribution       = 'Gauss';                           % Choose the expected type of noise distribution ('Gauss' or 'Rice')
settings.parameters.denoise.sigma              = 0;                                 % Specify the sigma of noise. If '0' noise will be estimated automatically
settings.parameters.denoise.verbose            = 0;

% Parameters for Wavelet denoising
settings.parameters.denoise.shrinkcheck        = false;
settings.parameters.denoise.NwLevel            = 4;
settings.parameters.denoise.wname              = 'db5';

% Parameters for Coupe's filters
settings.parameters.denoise.type               = 2;                                 % Choose which filter type to use of Coupé et al. algorithms (1 = AONLM; 2 = MR-ONLM; 3 = ONLM; 4 = ODCT; 5 = PRINLM)
settings.parameters.denoise.beta               = 0.2;
settings.parameters.denoise.patchradius        = 1;
settings.parameters.denoise.searchradius       = 3;

% Paramters for BM4D
settings.parameters.denoise.profile            = 'lc';                              % Choose profile for BM4D algorithm. (predefined: 'lc', 'np', 'mp' or 'manual')
settings.parameters.denoise.do_wiener          = true;
settings.parameters.denoise.N1                 = 8;      % (default: 4) cube has size (N1 x N1 x N3), needs to be of size x^2 
settings.parameters.denoise.Nstep              = 5;      % (default: 3) sliding step to process every next reference cube
settings.parameters.denoise.N2                 = 32;     % (default: 16) maximum number of similar cubes
settings.parameters.denoise.Ns                 = 11;     % (default: 11) search neighborhood size NsxNsxNs
settings.parameters.denoise.tau_match          = 0.1;    % (default: 0.1) threshold for the cube-distance (should be in [0,1])
settings.parameters.denoise.lambda_thr4D       = 1;      % (default: 2.7) hard-thresholding parameter, lower preserves more structures, higher smooths more
settings.parameters.denoise.N2_wiener          = 32;     % (default: 32) maximum number of similar cubes
settings.parameters.denoise.tau_match_wiener   = 0.05;   % (default: 0.05)
