function AcqWeightingFactor=getAcqWeightingFactors(settings, Lines, Partitions, amplitude)
% Function to get the acquisition weighting factors
% 
% ************************************************************************
% Version 1.0                                                  15.02.2019
% Falk Luesebrink              falk dot luesebrink at med dot ovgu dot de
% ************************************************************************

% compute 2D mask
w_lin = window(@hann, settings.parameters.nPE);
w_par = window(@hann, settings.parameters.nSlc);
[mask_row, mask_col] = meshgrid(w_par, w_lin);

% multiply by amplitude
w_mask = (mask_row .* mask_col * amplitude) ;

% compute AverageFactors
AcqWeightingFactor = w_mask(Lines, Partitions);
