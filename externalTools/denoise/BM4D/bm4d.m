function [y_est, sigma_est] = bm4d(z, settings)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  BM4D is an algorithm for attenuation of Gaussian or Rician noise in 
%  volumetric data. This algorithm reproduces the results from the articles:
%
%  [1] M. Maggioni, V. Katkovnik, K. Egiazarian, A. Foi, "A Nonlocal
%      Transform-Domain Filter for Volumetric Data Denoising and
%      Reconstruction", IEEE Trans. Image Process., vol. 22, no. 1, 
%      pp. 119-133, January 2013.  doi:10.1109/TIP.2012.2210725
%
%  [2] M. Maggioni, A. Foi, "Nonlocal Transform-Domain Denoising of 
%      Volumetric Data With Groupwise Adaptive Variance Estimation", 
%      Proc. SPIE Electronic Imaging 2012, San Francisco, CA, USA, Jan. 2012.
%
%
%  FUNCTION INTERFACE:
% 
%  INPUTS:
%     1) z               (3D array) : Noisy volume having any intensity range.
%     2) distribution        (char) : 'Gauss' --> z has Gaussian distribution
%                                   : 'Rice'  --> z has Rician distribution
%     3) sigma             (double) : Noise standard deviation, if unknown 
%                                     set it to 0 to enable noise estimation
%                                     (default is 0)
%     4) profile             (char) : 'lc' --> low complexity profile 
%                                   : 'np' --> normal profile 
%                                   : 'mp' --> modified profile
%                                     (default is 'mp')
%     5) do_wiener        (logical) : Perform collaborative Wiener filtering
%                                     (default is 1)
%     6) verbose          (logical) : 0 --> do not print output information
%                                     1 --> print information to screen
%                                     (default is 0)
%
%   Only z and distribution are required, all other inputs are optional. 
%   Optional inputs can be omitted, and assume their default value when 
%   set to 'empty' [].
%
%  OUTPUTS:
%     1) y_est           (3D array) : Final estimate (in the original range of z)
%     2) sigma_est       (3D array) : Voxel-wise standard deviation estimate
%
%
%  TYPICAL USAGE EXAMPLES:
%
%  Case: Noise with Gaussian distribution and known standard deviation
%   y_est = bm4d(z, 'Gauss', sigma);
%
%  Case: Noise with Rician distribution and unknown standard deviation
%   y_est = bm4d(z, 'Rice');
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Copyright (c) 2010-2015 Tampere University of Technology.
% All rights reserved.
% This work should only be used for nonprofit purposes.
%
% AUTHOR:
%     Matteo Maggioni, email: matteo.maggioni _at_ tut.fi
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Set parameters according to settings of processing pipeline (by Falk Luesebrink)
distribution = settings.parameters.denoise.distribution;
sigma        = settings.parameters.denoise.sigma;
profile      = settings.parameters.denoise.profile;
verbose      = settings.parameters.denoise.verbose;
do_wiener    = settings.parameters.denoise.do_wiener;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%  Arguments
if ~exist('z','var') || isempty(z)
    error('Argument "z" is required but was omitted.');
end
% Set Gaussian distribution by default
if ~exist('distribution','var') || isempty(distribution)
    distribution = 'Gauss';
    warning('Argument "distribution" is required but was omitted. Assuming Gaussian noise.');
end
% Enable adaptive noise estimation by default
if ~exist('sigma','var') || isempty(sigma) || (length(sigma)==1 && sigma==0)
    sigma = zeros(1,size(z,4));
end
% Enable modified profile by default
if ~exist('profile','var') || isempty(profile)
    profile = 'mp';
end
% Enable Wiener Filtering by default
if ~exist('do_wiener','var') || isempty(do_wiener)
    do_wiener = 1; 
end
% Enable verbose mode by default
if ~exist('verbose','var') || isempty(verbose)
    verbose = 1;
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Input check
%%%%
if ~isa(z,'numeric')
    error('Input "z" must be numeric.');
end
if ~isa(distribution,'char')
    error('Input "distribution" must be a string.');
end
if ~isa(sigma,'numeric')
    error('Input "sigma" must be numeric.');
end
if ~isa(profile,'char')
    error('Input "profile" must be a string.');
end
if ~(do_wiener==0 || do_wiener==1)
    error('Input "do_wiener" must be logical.');
end
if ~(verbose==0 || verbose==1)
    error('Input "verbose" must be logical.');
end
if ndims(z)>4
    error('Invalid input image dimension.');
end
if ~strcmpi(distribution,'Gauss') && ~strcmpi(distribution,'Rice')
    warning(['Invalid distribution argument "',distribution,'". Assuming Gaussian noise.']);
    distribution = 'Gauss';
end
if length(sigma)~=size(z,4)
    error('Input "sigma" must contain one value for each channel of "z".');
end
if ~strcmpi(profile,'lc') && ~strcmpi(profile,'np') && ~strcmpi(profile,'mp') && ~strcmpi(profile,'manual')
    warning(['Invalid profile argument "',profile,'". Assuming modified profile.']);
    profile = 'mp';
end
sigma = double(sigma);
z     = double(z);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Algorithm parameters
%%%%

% Transforms ('dct', 'dst', 'hadamard', or anything that is listed by 'help wfilters'):
transform_2D_HT_name       = 'bior1.5'; % 2-D spatial transform (Hard thresholding)
transform_3rd_dim_HT_name  = 'bior1.5'; % 1-D spatial transform
transform_4th_dim_HT_name  = 'haar';    % 1-D nonlocal transform
transform_2D_Wie_name      = 'dct';     % 2-D spatial transform (Wiener filtering)
transform_3rd_dim_Wie_name = 'dct';     % 1-D spatial transform
transform_4th_dim_Wie_name = 'haar';    % 1-D nonlocal transform


% Hard-thresholding (HT) parameters:
N1                   = 8;      % (default: 4) cube has size (N1 x N1 x N3), needs to be of size x^2 
N3                   = min(N1,size(z,3));   % cube has size (N1 x N1 x N3)
Nstep                = 6;      % (default: 3), sliding step to process every next reference cube
N2                   = 16;     % (default: 16) maximum number of similar cubes
Ns                   = 11;     % (default: 11) search neighborhood size NsxNsxNs
tau_match            = 0.1;    % (default: 0.1) threshold for the cube-distance (should be in [0,1])
lambda_thr4D         = 2.7;    % (default: 2.7) hard-thresholding parameter, lower preserves more structures, higher smooths more
thresholding         = 'hard'; % tresholding type
                               %  'hard' --> hard thresholding
                               %  'soft' --> soft thresholding
sharpening           = 1.0;    % sharpening alpha parameter (1 -> no sharpening)

% Wiener filtering parameters:
N1_wiener            = N1;
N3_wiener            = min(N1_wiener,size(z,3));
Nstep_wiener         = Nstep;
N2_wiener            = 32;
Ns_wiener            = Ns;
tau_match_wiener     = 0.05;     % (default: 0.05)

% Cube-matching parameters:
synchronous          = 0;  % if 1, the grouped cubes belong to the same slice
decLevel             = 0;  % decimation levels of the dyadic 2D wavelet transform

% Manually set parameters (by Falk Luesebrink)
if strcmp(profile,'manual')
    N1                   = settings.parameters.denoise.N1;                  % (default: 4) cube has size (N1 x N1 x N3), needs to be of size x^2 
    N3                   = min(settings.parameters.denoise.N1,size(z,3));   % cube has size (N1 x N1 x N3)
    Nstep                = settings.parameters.denoise.Nstep;               % (default: 3), sliding step to process every next reference cube
    N2                   = settings.parameters.denoise.N2;                  % (default: 16) maximum number of similar cubes
    Ns                   = settings.parameters.denoise.Ns;                  % (default: 11) search neighborhood size NsxNsxNs
    tau_match            = settings.parameters.denoise.tau_match;           % (default: 0.1) threshold for the cube-distance (should be in [0,1])
    lambda_thr4D         = settings.parameters.denoise.lambda_thr4D;        % (default: 2.7) hard-thresholding parameter, lower preserves more structures, higher smooths more

    % Wiener filtering parameters:
    N1_wiener            = settings.parameters.denoise.N1;
    N3_wiener            = min(N1_wiener,size(z,3));
    Nstep_wiener         = settings.parameters.denoise.Nstep;
    N2_wiener            = settings.parameters.denoise.N2_wiener;
    Ns_wiener            = settings.parameters.denoise.Ns;
    tau_match_wiener     = settings.parameters.denoise.tau_match_wiener;    % (default: 0.05)
end

% Modified profile parameters
if strcmpi(profile,'mp')
    N2               = 32;
    N1_wiener        = 5;
    N3_wiener        = min(N1_wiener,size(z,3));
    tau_match        = 0.2;
    tau_match_wiener = 0.1;
end

% Low complexity profile parameters
if strcmpi(profile,'lc')
    N2_wiener        = 16;
    Ns               = 7;
    Ns_wiener        = Ns;
    tau_match        = 0.04;         % default 0.04
    tau_match_wiener = 0.02;         % default 0.02
end

% Transform matrices
[Tfor, Tinv]     = getTransfMatrix(N1, transform_2D_HT_name, decLevel);
[Tfor3, Tinv3]   = getTransfMatrix(N3, transform_3rd_dim_HT_name, decLevel);
[TforW, TinvW]   = getTransfMatrix(N1_wiener, transform_2D_Wie_name, 0);
[Tfor3W, Tinv3W] = getTransfMatrix(N3_wiener, transform_3rd_dim_Wie_name, 0);

Tfor4 = cell(1,max(N2,N2_wiener));
Tinv4 = cell(1,max(N2,N2_wiener));
for hpow = 0:ceil(log2(max(N2,N2_wiener))),
    h = 2^hpow;
    [Tfor4rd, Tinv4rd] = getTransfMatrix(h, transform_4th_dim_HT_name, 0);
    Tfor4{h} = double(Tfor4rd);
    Tinv4{h} = double(Tinv4rd');
end

Tfor4W = cell(1,max(N2,N2_wiener));
Tinv4W = cell(1,max(N2,N2_wiener));
for hpow = 0:ceil(log2(max(N2,N2_wiener))),
    h = 2^hpow;
    [Tfor4rd, Tinv4rd] = getTransfMatrix(h, transform_4th_dim_Wie_name, 0);
    Tfor4W{h} = double(Tfor4rd);
    Tinv4W{h} = double(Tinv4rd');
end

% Power spectral densities of the noise (disabled in this script)
rootPSD = zeros(N1,N1,N3);
rootPSDW = zeros(N1_wiener,N1_wiener,N3_wiener);

% Maltab2MEX parameters
is_rician = strcmpi(distribution,'Rice');
hard = strcmpi(thresholding,'hard');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Filtering
%%%%

if verbose
    if ismatrix(z)
        fprintf('Filtering image of resolution %dx%dpx \n', size(z,1), size(z,2));
    elseif ndims(z)==3
        fprintf('Filtering volume of resolution %dx%dx%dpx \n', size(z,1), size(z,2), size(z,3));
    elseif ndims(z)==4
        fprintf('Filtering multichannel volume of resolution %dx%dx%dx%dpx \n', size(z,1), size(z,2), size(z,3), size(z,4));
    end
    if is_rician==1
        fprintf('\tRician noise distribution \n');
    else
        fprintf('\tGaussian noise distribution \n');
    end
    if any(sigma==0)
        fprintf('\tAdaptive noise estimation enabled \n');
    else
        fprintf('\tStandard deviation: %s \n', num2str(sigma));
    end
    if do_wiener==1
        fprintf('\tWiener filtering enabled \n');
    else
        fprintf('\tWiener filtering disabled \n');
    end
    fprintf(['\tParameter profile "',profile,'" \n']);
end

% Hard thresholding
basic = tic;
[y_hat, sigma_hat] = bm4d_thr_mex(z, Nstep, N1, N2, N3,...
    lambda_thr4D, tau_match, (Ns-1)/2, synchronous, ...
    sigma, sharpening, hard, is_rician, do_wiener, Tfor, Tinv', Tfor3, ...
    Tinv3', Tfor4, Tinv4, rootPSD, transform_2D_HT_name, ...
    transform_3rd_dim_HT_name, transform_4th_dim_HT_name );
if verbose
    fprintf('Basic estimate completed (%.1fs) \n', toc(basic));
end

% Wiener filtering
if do_wiener
    wiener = tic;
    [y_est, sigma_est] = bm4d_wie_mex(z, y_hat, Nstep_wiener, N1_wiener, ...
        N2_wiener, N3_wiener, tau_match_wiener, (Ns_wiener-1)/2, ...
        synchronous, sigma, is_rician, TforW, TinvW', Tfor3W, Tinv3W', Tfor4, Tinv4, rootPSDW, ...
        transform_2D_Wie_name, transform_3rd_dim_Wie_name, transform_4th_dim_Wie_name );
    if verbose
        fprintf('Final estimate completed (%.1fs) \n', toc(wiener))
    end
else
    y_est     = y_hat;
    sigma_est = sigma_hat;
end

if verbose
    fprintf('Total execution time: %.1f sec \n', toc(basic));
end

return





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Some auxiliary functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




function [Transf_Matrix, invTransf_Matrix] = getTransfMatrix (N, transform_type, Nden)
%
% Create forward and inverse transform matrices, which allow for perfect
% reconstruction. The forward transform matrix is normalized so that the
% l2-norm of each basis element is 1.
%
% [Tforward, Tinverse] = getTransfMatrix (N, transform_type, dec_levels)
%
%  INPUTS:
%
%   N           --> Size of the transform (for wavelets, must be 2^K)
%
%   transform_type  --> 'dct', 'dst', 'hadamard', or anything that is
%                       listed by 'help wfilters' (bi-orthogonal wavelets)
%                       'DCrand' -- an orthonormal transform with a DC and all
%                       the other basis elements of random nature
%
%   dec_levels      --> If a wavelet transform is generated, this is the
%                       desired decomposition level. Must be in the
%                       range [0, log2(N)-1], where "0" implies
%                       full decomposition.
%
%  OUTPUTS:
%
%   Tforward        --> (N x N) Forward transform matrix
%
%   Tforward        --> (N x N) Inverse transform matrix
%

if ~exist('Nden','var')
    Nden = 0;
end

if N == 1,
    Transf_Matrix = 1;
elseif strcmp(transform_type, 'dct') == 1,
    if N<=16
        load Transforms.mat TDCTFor
        Transf_Matrix    = TDCTFor{N};
    else
        Transf_Matrix    = dct(eye(N));
    end
elseif strcmp(transform_type, 'dst') == 1,
    if N<=16
        load Transforms.mat TDSTFor
        Transf_Matrix    = TDSTFor{N};
    else
        Transf_Matrix    = dst(eye(N));
    end
elseif strcmp(transform_type, 'DCrand') == 1,
    x = randn(N); x(1:end,1) = 1; [Q,R] = qr(x);
    if (Q(1) < 0),
        Q = -Q;
    end;
    Transf_Matrix = Q';
elseif strcmp(transform_type, 'hadamard') == 1,
    Transf_Matrix    = hadamard(N);
else %% wavelet transform

    if (strcmpi(transform_type,'bior1.5') || strcmpi(transform_type,'haar')) && N<=64
        if strcmpi(transform_type,'bior1.5')
            load Transforms.mat TBiorFor
            TR = TBiorFor;
        else
            load Transforms.mat THaarFor
            TR = THaarFor;
        end
        Transf_Matrix = TR{N};
    else
        %%% Set periodic boundary conditions, to preserve bi-orthogonality
        dwtmode('per','nodisp');

        [LO_D, HI_D, LO_R, HI_R] = wfilters(transform_type);
        for i = 1:N
            Transf_Matrix(i,:)=waverec(circshift([1 zeros(1,N-1)],[Nden i-1]), ...
                2.^[Nden Nden:log2(N)], LO_D, -HI_D);  %% construct transform matrix
        end
    end
    
end

%%% Normalize the basis elements
Transf_Matrix = (Transf_Matrix' * diag(sqrt(1./sum(Transf_Matrix.^2,2))))';

%%% Compute the inverse transform matrix
invTransf_Matrix = inv(Transf_Matrix);

return;
