function parsave(pathandfilename, data)
% Function to save data using parfor
%
% ************************************************************************
% Version 1.0                                                  15.02.2019
% Falk Luesebrink              falk dot luesebrink at med dot ovgu dot de
% ************************************************************************
    
save(pathandfilename, 'data', '-v7.3', '-nocompression');