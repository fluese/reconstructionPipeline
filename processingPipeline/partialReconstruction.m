function kspace=partialReconstruction(settings, kspace)
% This function synthesizes data in case of k space undersampling with
% partial Fourier.
%
% Possible methods:
% 1) Zero filling
% 2) Conjugate synthesize
% 3) POCS
% 4) Homodyne
%
%
% ************************************************************************
% Version 1.0                                                  15.02.2019
% Falk Luesebrink              falk dot luesebrink at med dot ovgu dot de
% ************************************************************************

if settings.options.kspace.partialFourier
    if strcmp(settings.parameters.kspace.partialFourier, 'zero')
        full_kspace = complex(zeros(settings.parameters.nRO, settings.parameters.nPE, settings.parameters.nSlc,'single'));
        full_kspace(settings.parameters.uRO+1:end,settings.parameters.uPE+1:end,settings.parameters.uSlc+1:end) = kspace;
    elseif strcmp(settings.parameters.kspace.partialFourier, 'conj')
%         fprintf('| Conjugate synthesize currently not implemented.\n')
    elseif strcmp(settings.parameters.kspace.partialFourier, 'pocs')
%         fprintf('| POCS currently not implemented.\n')
    elseif strcmp(settings.parameters.kspace.partialFourier, 'homo')
%         fprintf('| Homodyne currently not implemented.\n')
    elseif uRO == 0 && uPE == 0 && uSlc == 0
%         fprintf('| No partial Fourier reconstruction needed. Proceeding without it.\n')
    else
%         fprintf('| Invalid method for partial Fourier reconstruction specified. Proceeding without filling k space.\n')
        return;
    end
    kspace=full_kspace;
else
    return;
end