function settings=statusData(name, settings, type)
% Function to check whether data have already been processed.
%
% ************************************************************************
% Version 1.0                                                  15.02.2019
% Falk Luesebrink              falk dot luesebrink at med dot ovgu dot de
% ************************************************************************

if strcmp(type, 'kspace')
    if settings.options.force || settings.options.kspace.force
        fprintf('| **** WARNING! Forcing to write data.\n')
        fprintf('| \n');
        settings.parameters.kspace.processed = false;
        settings.options.xspace.force        = true;
        settings.options.imspace.force       = true;
        settings.options.NIfTI.force         = true; 
    elseif exist([name.namePath '/' name.kspace '/Slice_1/Channel_1.mat'], 'file')
        settings.parameters.kspace.processed = true;
    else
        settings.parameters.kspace.processed = false;
    end

    if ~exist([name.namePath '/' name.kspace '/'], 'dir')
        mkdir([name.namePath '/' name.kspace '/']);
    end
    
    if ~exist([name.namePath '/' name.kspace '/Slice_1/'], 'dir')
        for Slc = 1:settings.parameters.pSlc
            mkdir([name.namePath '/' name.kspace '/Slice_' num2str(Slc) '/']);
        end
    end
elseif strcmp(type, 'xspace')
    if settings.options.force || settings.options.xspace.force
        fprintf('| **** WARNING! Forcing to write data.\n')
        fprintf('| \n');
        settings.parameters.xspace.processed = false;
        settings.options.imspace.force       = true;
        settings.options.NIfTI.force         = true;
    elseif exist([name.namePath '/' name.xspace '/Slice_1.mat'], 'file')
        settings.parameters.xspace.processed = true;
    else
        settings.parameters.xspace.processed = false;
    end

    if ~exist([name.namePath '/' name.xspace '/'], 'dir')
        mkdir([name.namePath '/' name.xspace '/']);
    end
    
    if ~exist([name.namePath '/' name.xspace '/Slice_1/'], 'dir')
        for Slc = 1:settings.parameters.nSlc
            mkdir([name.namePath '/' name.xspace '/Slice_' num2str(Slc) '/']);
        end
    end
elseif strcmp(type, 'imspace')
    if settings.options.force || settings.options.imspace.force
        fprintf('| **** WARNING! Forcing to write data.\n')
        fprintf('| \n');
        settings.parameters.imspace.processed = false;
        settings.options.NIfTI.force          = true;
    elseif exist([name.namePath '/' name.imspace '/Slice_1.mat'], 'file')
        settings.parameters.imspace.processed = true;
    else
        settings.parameters.imspace.processed = false;
    end

    if ~exist([name.namePath '/' name.imspace '/'], 'dir')
        mkdir([name.namePath '/' name.imspace '/']);
    end
elseif strcmp(type, 'NIfTI')
    if settings.options.force || settings.options.NIfTI.force
        fprintf('| **** WARNING! Forcing to write data.\n')
        fprintf('| \n');
        settings.parameters.NIfTI.processed = false;
    elseif exist([name.namePath '/' name.nameFile '.nii.gz'], 'file') || exist([name.namePath '/' name.nameFile '.nii'], 'file')
        settings.parameters.NIfTI.processed = true;
    else
        settings.parameters.NIfTI.processed = false;
    end
end