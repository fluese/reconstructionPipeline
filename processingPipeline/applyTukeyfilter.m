function kspace=applyTukeyfilter(settings, kspace)
% Function to filter k space data with Tukey window to reduce Gibb's
% ringing.
%
% ************************************************************************
% Version 1.0                                                  15.02.2019
% Falk Luesebrink              falk dot luesebrink at med dot ovgu dot de
% ************************************************************************

if settings.options.kspace.tukeyWindow
    filterWindow = window3_tukeytaper([settings.parameters.pRO settings.parameters.pPE settings.parameters.pSlc], @tukeywin, settings.parameters.kspace.tukey);
    kspace       = kspace .* filterWindow;
end