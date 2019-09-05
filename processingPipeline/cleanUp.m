function cleanUp(name, settings, numFiles)
% Function to clean up after processing.
%
% ************************************************************************
% Version 1.0                                                  15.02.2019
% Falk Luesebrink              falk dot luesebrink at med dot ovgu dot de
% ************************************************************************
fprintf('\n');
fprintf('===============================================================\n');
fprintf('| Cleaning up.\n');
fprintf('===============================================================\n');

settings.endTime = datestr(datetime('now'));

fprintf('| Writing settings file in subject folder.\n');
save([name.namePath '/Settings_' name.nameFile '.mat'], 'settings', '-v7.3')

if ~settings.options.keepFiles && ~settings.options.saveIntermediate
    fprintf('| Removing all intermediately written files.\n');
    if exist([name.namePath '/' name.kspace], 'dir')
        rmdir([name.namePath '/' name.kspace], 's');
    end
    if exist([name.namePath '/' name.xspace], 'dir')
        rmdir([name.namePath '/' name.xspace], 's');
    end
    if exist([name.namePath '/' name.imspace],'dir')
        rmdir([name.namePath '/' name.imspace],'s');
    end
elseif settings.options.keepFiles
    fprintf('| Keeping all intermediately written files.\n');
end

if numFiles == length(name.FilenameS)
    poolobj = gcp('nocreate');
    if ~isempty(poolobj)
        fprintf('| ');
        delete(gcp('nocreate'));
    end
end

fprintf('===============================================================\n');
fprintf('| Reconstruction finished at %s\n', datestr(datetime('now')));
fprintf('===============================================================\n\n');