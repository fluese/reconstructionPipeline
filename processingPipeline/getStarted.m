function name=getStarted
% Function to get the processing pipeline started.
% 
% ************************************************************************
% Version 1.0                                                  15.02.2019
% Falk Luesebrink              falk dot luesebrink at med dot ovgu dot de
% ************************************************************************

[name.FilenameS, name.Filepath] = uigetfile('*.dat', 'Select raw data files to be reconstructed.','MultiSelect', 'on');

% Manual setup in case of -nodisplay!
% name.Filepath  	    = '/pool/falk/rawdata/250um/';
% name.FilenameS      = cell(8,1);
% name.FilenameS{1,1} = 'meas_MID168_mprage_xpace_0_25iso_testing_FID11814.dat';
% name.FilenameS{2,1} = 'meas_MID74_mprage_xpace_0_25iso_testing_FID14411.dat';
% name.FilenameS{3,1} = 'meas_MID75_mprage_xpace_0.25iso_testing_FID14412.dat';
% name.FilenameS{4,1} = 'meas_MID40_mprage_xpace_0_25iso_testing_FID16596.dat';
% name.FilenameS{5,1} = 'meas_MID24_mprage_xpace_0_25iso_testing_FID19885.dat';
% name.FilenameS{6,1} = 'meas_MID25_mprage_xpace_0_25iso_testing_FID19886.dat';
% name.FilenameS{7,1} = 'meas_MID31_mprage_xpace_0.25iso_testing_FID24165.dat';
% name.FilenameS{8,1} = 'meas_MID32_mprage_xpace_0.25iso_testing_FID24166.dat';
if ~iscell(name.FilenameS)           % If one file was specified only, change it to a cell, in order to make it work.
    tmp = name.FilenameS;
    name.FilenameS=cell(1,1);
    name.FilenameS{1,1} = tmp;
    clear tmp
end

if name.FilenameS{1} == 0
    error('No files selected.')
end
