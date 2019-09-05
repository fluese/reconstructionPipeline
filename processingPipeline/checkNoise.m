function isNoise=checkNoise(preset)
% Function to check whether noise data has been processed.
%
% ************************************************************************
% Version 1.0                                                  15.02.2019
% Falk Luesebrink              falk dot luesebrink at med dot ovgu dot de
% ************************************************************************

if preset.usePreset && strcmp(preset.file, 'noise')
    isNoise = true;
    
    poolobj = gcp('nocreate');
    if ~isempty(poolobj)
        fprintf('| ');
        delete(gcp('nocreate'));
    end
    
    fprintf('===============================================================\n');
    fprintf('| Reconstruction finished at %s\n', datestr(datetime('now')));
    fprintf('===============================================================\n');
else
    isNoise = false;
end