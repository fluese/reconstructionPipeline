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
fpath = uigetdir(userdir,'Select folder where data shall be saved into.');    
spath = uigetdir(userdir,'Select folder where spm is located at.');

% Setup manually in case of -nodisplay
% fpath = '/pool/falk/reconstruction';
% spath = '/home/luesebri/software/spm12';

A = regexp(fileread('./processingPipeline/setupPipeline.m'), '\n', 'split');
if ispc
    tmpLocation = ['settings.parameters.path      = ''' fpath '\'';                  % Set path where data shall be written in.'];
    A{73} = sprintf('%s', tmpLocation);
    if spath == 0
        tmpLocation = 'settings.parameters.post.path = ''none'';               % Set path to tissue probability model';
    else
        tmpLocation = ['settings.parameters.post.path = ''' spath '\tpm\TPM.nii'';       % Set path to tissue probability model'];
    end
    A{74} = sprintf('%s', tmpLocation);
else
    tmpLocation = ['settings.parameters.path      = ''' fpath '/'';                  % Set path where data shall be written in.'];
    A{73} = sprintf('%s', tmpLocation);
    if spath == 0
        tmpLocation = 'settings.parameters.post.path = ''none'';               % Set path to tissue probability model';
    else
        tmpLocation = ['settings.parameters.post.path = ''' spath '/tpm/TPM.nii'';       % Set path to tissue probability model'];
    end
    A{74} = sprintf('%s', tmpLocation);
end
fid = fopen('./processingPipeline/setupPipeline.m', 'w');
fprintf(fid, '%s\n', A{:});
fclose(fid);
