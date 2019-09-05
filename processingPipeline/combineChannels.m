function imspace=combineChannels(name, settings, Chunk, xspace)
% Function to combine channels by (root) sum of squares or adaptive combine.
%
% ************************************************************************
% Version 1.0                                                  15.02.2019
% Falk Luesebrink              falk dot luesebrink at med dot ovgu dot de
% ************************************************************************

Index = ((Chunk-1)*settings.parameters.stepSize)+1 : Chunk*settings.parameters.stepSize;

% If index requests a slice beyond the maximum number of slices, remove that element
if Index(end) > settings.parameters.nSlc
    entriesToRemove = Index > settings.parameters.nSlc;
    Index(entriesToRemove) = [];
    removedEntries = sum(entriesToRemove);
else
    removedEntries = 0;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Channel combination
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if strcmp(settings.options.xspace.combineChannels, 'adaptComb')
    % Adaptive combine
    if ~settings.options.procMemory
        xspace=loadData(name, settings, Index, Chunk, 'adaptComb');
    end
    imspace = adaptiveCombine(settings, xspace);
elseif strcmp(settings.options.xspace.combineChannels, 'SoS')
    % Root Sum of squares
    if ~settings.options.procMemory
        xspace=loadData(name, settings, Index, Chunk, 'SoS', removedEntries);
    end
    xspace  = permute(xspace,[2 3 4 1]);
    imspace = squeeze(sum(abs(xspace.^2),4)).^(1/2);
end

if settings.options.saveIntermediate && settings.options.procMemory
    saveData(name, settings, imspace, Index, Chunk, 'imspace', removedEntries)
elseif ~settings.options.procMemory
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Save data
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    saveData(name, settings, imspace, Index, Chunk, 'imspace', removedEntries)
    imspace=0;
end