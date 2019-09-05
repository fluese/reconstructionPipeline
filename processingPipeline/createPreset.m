function createPreset(name)
% Function to create a preset of settings specified in 'setOptions.m' and
% 'setParamters.m'.
%
% ************************************************************************
% Version 1.0                                                  15.02.2019
% Falk Luesebrink              falk dot luesebrink at med dot ovgu dot de
% ************************************************************************
settings = setOptions;
settings = setParameters(settings);

Path = fileparts(mfilename('fullpath'));
save([Path '/presets/' name '.mat'], 'settings', '-v7.3');