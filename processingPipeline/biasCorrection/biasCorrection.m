function biasCorrection(settings, filename)
% This is a script to run SPM's bias field correction without the GUI and
% to circumvent limitations made by it (e.g. lowest option of FWHM is 30).
%
% Falk Lüsebrink (22.06.2016)
% 
% Minor Update
% Falk Lüsebrink (22.11.2017)

% Set path to SPM's tissue probability model file.
TPMpath = settings.parameters.post.path;

% Here you can set FWHM and regularization for SPM's bias field correction.
%
% What works best will depend on the amount of bias artifact in your data.
% If there is no bias artifact, you would use very heavy regularisation
% because there should be nothing to correct. If your data are heavily
% corrupted by bias artifact, then you'd use much less regularisation in
% order that the model has more flexibility.
% The FWHM option also depends on the nature of the artifact. If it is very
% low frequency, you'd be better off with a broad FWHM. If it contains high
% frequency, you may need a smaller FWHM.
fwhm = settings.parameters.post.fwhm; % Default: 60
reg  = settings.parameters.post.reg;  % Default: 0.001

% This is the option to save a bias corrected version of your images and/or
% the estimated bias field. [0,0] disables both, [0,1] enables writing of
% corrected version of your image, [1,0] enables writing of estilmated bias
% field and [1,1] enables both. Default: [1,1]
write = [0,1];

% Set sampling distance between points in mm. Default: 3.
samp = settings.parameters.post.samp;

jobs = cell(1,1);

% Load batchfile, extract job and clear batchfile afterwards
load defjob.mat
jobs{1,1} = matlabbatch{1,1}.spm.spatial.preproc; clear matlabbatch;

% Set path to SPM model and disable writing of segmentation files,
% which should reduce processing time. [Unconfirmed]
% In this case segmentation files of GM, WM, and CSF are written, as they
% are used to generate a brain mask.
if settings.options.post.segmentation
    for j=1:3
    jobs{1,1}.tissue(1,j).tpm = {TPMpath};
    jobs{1,1}.tissue(1,j).native = [1,0];
    end
else
    for j=1:3
    jobs{1,1}.tissue(1,j).tpm = {TPMpath};
    jobs{1,1}.tissue(1,j).native = [0,0];
    end
end

for j=4:6
jobs{1,1}.tissue(1,j).tpm = {TPMpath};
jobs{1,1}.tissue(1,j).native = [0,0];
end


% Change parameters of job specified above
jobs{1,1}.channel.biasfwhm = fwhm;
jobs{1,1}.channel.biasreg  = reg;
jobs{1,1}.channel.write    = write;
jobs{1,1}.warp.samp        = samp;

% Setup path and volume to be processed
jobs{1,1}.channel.vols = {filename};

% Run SPM bias field correction with job as settings.parameters. Creates a new
% nifti with '_corrected' as appedix in the same folder as input data.
spm_preproc_run_fl(jobs{1,1},'run');