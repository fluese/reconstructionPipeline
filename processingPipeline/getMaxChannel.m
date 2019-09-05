function settings=getMaxChannel(name, settings, xspace)
% Function to estimate the channel with highest intesity for adaptive
% combine.
%
% ************************************************************************
% Version 1.0                                                  15.02.2019
% Falk Luesebrink              falk dot luesebrink at med dot ovgu dot de
% ************************************************************************

if settings.options.xspace.getMaxChannel && strcmp(settings.options.xspace.combineChannels, 'adaptComb')
    fprintf('| The channel with the highest intensity will be retrieved for phase correction.\n');
    if settings.options.procMemory
        [~,settings.parameters.maxchannel]=max(sum(abs(xspace(:,:)),2));
    elseif settings.options.procParallel
        fprintf('| The data will be split in %i chunks and %i chunks are processed in parallel.\n', settings.parameters.nChunks, settings.parameters.combineWorkers);
        value = zeros(settings.parameters.nChunks,1);
        coil = zeros(settings.parameters.nChunks,1);
        stepSize=settings.parameters.stepSize;
        
        parfor Chunk = 1:settings.parameters.nChunks
            Index = ((Chunk-1)*stepSize)+1 : Chunk*stepSize;
            Im=loadData(name, settings, Index, Chunk, 'maxChannel');
            [value(Chunk),coil(Chunk)]=max(sum(abs(Im(:,:)),2));
        end
        clear Im
        [~,idx] = max(value);
        settings.parameters.maxchannel = coil(idx);
    elseif ~settings.options.procParallel
        Im=loadData(name, settings, [], [], 'maxChannel');
        [~,settings.parameters.maxchannel]=max(sum(abs(Im(:,:)),2));
        clear Im
    end
    fprintf('| The channel with the highest intensity is %i.\n', settings.parameters.maxchannel);
    fprintf('| \n');
else
    settings.parameters.maxchannel=1;
end