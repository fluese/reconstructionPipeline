function saveData(name, settings, data, Index, Chunk, type, removedEntries)
% Function to save input data to disk by type of being processed.
%
% ************************************************************************
% Version 1.0                                                  15.02.2019
% Falk Luesebrink              falk dot luesebrink at med dot ovgu dot de
% ************************************************************************

if strcmp(type, 'kspace') % Save data in k space.
    for Slc = 1: settings.parameters.nSlc
        path = [name.namePath '/' name.kspace '/Slice_' num2str(Slc) '/Channel_' num2str(Chunk) '.mat'];
        parsave(path, squeeze(data(:,:,Slc)))
    end
elseif strcmp(type, 'preproc') % Save data in k space after pre-processing data
    for Channel = 1:settings.parameters.nCh
%         if settings.options.saveIntermediate && settings.options.procParallel
%              for Slc = Index
%                 path = [name.namePath '/' name.kspace '/Slice_' num2str(Slc) '/Channel_' num2str(Channel) '.mat'];
%                 parsave(path, squeeze(data(Channel,:,:,Slc-(Chunk-1)*settings.parameters.stepSize)));
%             end
%         elseif settings.options.saveIntermediate
%             for Slc = 1:settings.parameters.pSlc
%                 path = [name.namePath '/' name.kspace '/Slice_' num2str(Slc) '/Channel_' num2str(Channel) '.mat'];
%                 parsave(path, squeeze(data(Channel,:,:,Slc)));
%             end
%         else
            for Slc = Index
                path = [name.namePath '/' name.kspace '/Slice_' num2str(Slc) '/Channel_' num2str(Channel) '.mat'];
%                 disp(['Trying to write: ' path]);
                parsave(path, squeeze(data(Channel,:,:,Slc-(Chunk-1)*settings.parameters.stepSize)));
            end
%         end
    end
elseif strcmp(type, 'xspace') % Save data in x space (and do inverse FFT).
    % Fourier transform k space into x space
    fft_dims           = [1 2 3];
    for f = fft_dims
        data = ifftshift(ifft(fftshift(data,f),[],f),f);
    end
    
    for Slc = 1:settings.parameters.nSlc
        path = [name.namePath '/' name.xspace '/Slice_' num2str(Slc) '/Channel_' num2str(Chunk) '.mat'];
        parsave(path, squeeze(data(:,:,Slc)));
    end
elseif strcmp(type, 'sepChannels') % Save data from separated channels into a single file per slice (hopefully reducing disk IO)
    if settings.options.procParallel
        nCh      = settings.parameters.nCh;
        nRO      = settings.parameters.nRO;
        nPE      = settings.parameters.nPE;
        nSlc     = settings.parameters.nSlc;
        namePath = name.namePath;
        nxspace  = name.xspace;
        
        parfor Slc = 1:nSlc
            Im=complex(zeros(nCh,nRO,nPE, 'single'));
            for Ch = 1:nCh
                kspace = load([namePath '/' nxspace '/Slice_' num2str(Slc) '/Channel_' num2str(Ch) '.mat']);
                Im(Ch, :, :) = kspace.data;
            end
            parsave([namePath '/' nxspace '/Slice_' num2str(Slc) '.mat'], Im);
            rmdir([namePath '/' nxspace '/Slice_' num2str(Slc)], 's');
        end
    else
        for Slc = 1:settings.parameters.nSlc
            Im=complex(zeros(settings.parameters.nCh,settings.parameters.nRO,settings.parameters.nPE, 'single'));
            for Ch = 1:settings.parameters.nCh
                kspace = load([name.namePath '/' name.xspace '/Slice_' num2str(Slc) '/Channel_' num2str(Ch) '.mat']);
                Im(Ch, :, :) = kspace.data;
            end
            parsave([name.namePath '/' name.xspace '/Slice_' num2str(Slc) '.mat'], Im);
            rmdir([name.namePath '/' name.xspace '/Slice_' num2str(Slc)], 's');
        end
    end
elseif strcmp(type, 'saveIntermediate')
    if settings.options.procParallel
        nCh      = settings.parameters.nCh;
        nRO      = settings.parameters.nRO;
        nPE      = settings.parameters.nPE;
        nSlc     = settings.parameters.nSlc;
        namePath = name.namePath;
        nxspace  = name.xspace;
        
        parfor Slc = 1:nSlc
            Im=complex(zeros(nCh,nRO,nPE, 'single'));
            for Ch = 1:nCh
                Im(Ch, :, :) = data(Ch,:,:,Slc);
            end
            parsave([namePath '/' nxspace '/Slice_' num2str(Slc) '.mat'], Im);
            rmdir([namePath '/' nxspace '/Slice_' num2str(Slc)], 's');
        end
    else
        for Slc = 1:settings.parameters.nSlc
            Im=complex(zeros(settings.parameters.nCh,settings.parameters.nRO,settings.parameters.nPE, 'single'));
            for Ch = 1:settings.parameters.nCh
                Im(Ch, :, :) = data(Ch,:,:,Slc);
            end
            parsave([name.namePath '/' name.xspace '/Slice_' num2str(Slc) '.mat'], Im);
            rmdir([name.namePath '/' name.xspace '/Slice_' num2str(Slc)], 's');
        end
    end    
elseif strcmp(type, 'imspace') % Save data in image space.
    path = [name.namePath '/' name.imspace '/'];
    
    if ~exist([name.namePath '/' name.imspace '/'], 'dir')
        mkdir([name.namePath '/' name.imspace '/']);
    end
    
    if strcmp(settings.options.xspace.combineChannels,'SoS') || strcmp(settings.options.xspace.combineChannels,'adaptComb') && strcmp(settings.parameters.xspace.combineType,'bySlice')
        for Slc = Index
            parsave([path 'Slice_' num2str(Slc) '.mat'], squeeze(data(:,:,Slc-((Chunk-1)*(length(Index)+removedEntries))))); % Slc-((nChunk-1)*settings.parameters.stepSize))) results in indicies from 1 to stepSize
        end
    elseif strcmp(settings.options.xspace.combineChannels,'adaptComb') && strcmp(settings.parameters.xspace.combineType,'byChunk')
        if Chunk == 1
            for Slc = Index
                parsave([path 'Slice_' num2str(Slc) '.mat'], data(:,:,Slc));
            end
        else
            for Slc = Index
                parsave([path 'Slice_' num2str(Slc) '.mat'], data(:,:,Slc-((Chunk-1)*settings.parameters.stepSize)));
            end
        end
    end
end