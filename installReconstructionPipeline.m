function installReconstructionPipeline
ipath = fileparts(mfilename('fullpath'));

% Add path to Matlab path
addpath(...
    [ipath '/']);

% Add path with subfolders to Matlab path
addpath(...
    genpath([ipath '/externalTools/']), ...
    genpath([ipath '/processingPipeline/']));

if ispc
  userdir = getenv('USERPROFILE');
else
  userdir = getenv('HOME');
end
fpath = uigetdir(userdir,'Select folder were data shall be saved into.');    
spath = uigetdir(userdir,'Select folder were spm is located at.');

% Setup manually in case of -nodisplay
% fpath = '/pool/falk/reconstruction';
% spath = '/home/luesebri/software/spm12';

A = regexp(fileread('./processingPipeline/setupPipeline.m'), '\n', 'split');
if ispc
    tmpLocation = ['settings.parameters.path      = ''' fpath '\'';                  % Set path were data shall be written in.'];
    A{68} = sprintf('%s', tmpLocation);
    if isempty(spath)
        % Put somewhere in setupPipeline that SPM is not installed!
    else
        tmpLocation = ['settings.parameters.post.path = ''' spath '\tpm\TPM.nii'';       % Set path to tissue probability model'];
        A{69} = sprintf('%s', tmpLocation);
    end
else
    tmpLocation = ['settings.parameters.path      = ''' fpath '/'';                  % Set path were data shall be written in.'];
    A{68} = sprintf('%s', tmpLocation);
    if isempty(spath)
        % Put somewhere in setupPipeline that SPM is not installed!
    else
        tmpLocation = ['settings.parameters.post.path = ''' spath '/tpm/TPM.nii'';       % Set path to tissue probability model'];
        A{69} = sprintf('%s', tmpLocation);
    end
end
fid = fopen('./processingPipeline/setupPipeline.m', 'w');
fprintf(fid, '%s\n', A{:});
fclose(fid);
