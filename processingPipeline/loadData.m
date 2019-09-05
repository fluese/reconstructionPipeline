function data=loadData(name, settings, Index, Chunk, type, removedEntries)
% Load data by type of processing.
%
% ************************************************************************
% Version 1.0                                                  15.02.2019
% Falk Luesebrink              falk dot luesebrink at med dot ovgu dot de
% ************************************************************************

nameVol  = name.namePath;
nSlc     = settings.parameters.nSlc;
nxspace  = name.xspace;
nkspace  = name.kspace;

if strcmp(type, 'adaptComb') % Loading data to combine channels with adaptive combine
    if strcmp(settings.parameters.xspace.combineType, 'bySlice') % If data shall be processed slice by slice.
        data = complex(zeros(settings.parameters.nCh, settings.parameters.nRO, settings.parameters.nPE, length(Index),'single'));
        for Slc = Index
            tmp=load([nameVol '/' nxspace '/Slice_' num2str(Slc) '.mat']);
            data(:,:,:,Slc-((Chunk-1)*settings.parameters.stepSize)) = tmp.data;
        end
    elseif strcmp(settings.parameters.xspace.combineType, 'byChunk') % If data shall be processed by chunks of slices.
        if Chunk == 1
            % Update Index due to offset, which is needed to avoid effects
            % from voxels at the border.
            % idx_offset = Index(end)+offset;
            % uIndex = [Index idx_offset];
            
            data = complex(zeros(settings.parameters.nCh, settings.parameters.nRO, settings.parameters.nPE, length(Index),'single'));
            for Slc = Index
                tmp=load([nameVol '/' nxspace '/Slice_' num2str(Slc) '.mat']);
                data(:,:,:,Slc-((Chunk-1)*length(Index))) = tmp.data;
            end
        else
            % Update Index due to offset, which is needed to avoid effects
            % from voxels at the border.
            % idx_offset = [Index(1)-settings.parameters.xspace.offset Index(end)+offset];
            % uIndex = idx_offset(1):idx_offset(2);
            
            % If index would request a slice beyond the maximum number of
            % slices, remove that element.
            if Index(end) > settings.parameters.nSlc     
                entriesToRemove = Index > settings.parameters.nSlc;
                Index(entriesToRemove) = [];
            end
            
            data = zeros(settings.parameters.nCh, settings.parameters.nRO, settings.parameters.nPE, length(Index), 'single');
            for Slc = Index
                tmp=load([nameVol '/' nxspace '/Slice_' num2str(Slc) '.mat']);
                data(:,:,:,Slc-((Chunk-1)*length(Index))) = tmp.data;
            end
        end
    end
elseif strcmp(type, 'SoS') % Loading data to combine channels with (root) sum of squares.
    data = complex(zeros(settings.parameters.nCh, settings.parameters.nRO, settings.parameters.nPE, length(Index), 'single'));
    for Slc = Index
        tmp=load([nameVol '/' nxspace '/Slice_' num2str(Slc) '.mat']);
        data(:,:,:,Slc-((Chunk-1)*(length(Index)+removedEntries))) = tmp.data;
    end
elseif strcmp(type, 'saveIntermediate') % Loading data to combine channels with (root) sum of squares.
    data=complex(zeros(settings.parameters.nCh, settings.parameters.nRO, settings.parameters.nPE, settings.parameters.nSlc, 'single'));
    for Slc = Index
        tmp=load([nameVol '/' nxspace '/Slice_' num2str(Slc) '.mat']);
        data(:,:,:,Slc) = tmp.data;
    end
elseif strcmp(type, 'maxChannel') % Loading data to retrieve the coil with the maximum intensity for phase correction in adaptive combine
    if settings.options.procParallel
        data = complex(zeros(settings.parameters.nCh, settings.parameters.nRO, settings.parameters.nPE, length(Index), 'single'));

        % If index requests a slice beyond the maximum number of slices, remove that element
        if Index(end) > settings.parameters.nSlc
            entriesToRemove = Index > settings.parameters.nSlc;
            Index(entriesToRemove) = [];
        end
        
        for Slc = Index
            tmp=load([nameVol '/' nxspace '/Slice_' num2str(Slc) '.mat']);
            data(:,:,:,Slc-((Chunk-1)*length(Index))) = tmp.data;
        end
    else
        data = complex(zeros(settings.parameters.nCh, settings.parameters.nRO, settings.parameters.nPE, settings.parameters.nSlc, 'single'));
        for k = 1:nSlc
            tmp=load([nameVol '/' nxspace '/Slice_' num2str(k) '.mat']);
            data(:,:,:,k) = tmp.data;
        end
    end
elseif strcmp(type, 'kspace') % Loading k space data to be reconstructed
    data = complex(zeros(settings.parameters.nCh, settings.parameters.nRO, settings.parameters.nPE, length(Index), 'single'));
    for Channel = 1:settings.parameters.nCh
        for Slc = Index
            tmp=load([nameVol '/' nkspace '/Slice_' num2str(Slc) '/Channel_' num2str(Channel) '.mat']);
            data(Channel,:,:,Slc-(Chunk-1)*length(Index)) = tmp.data;
        end
    end
elseif strcmp(type, 'kspace_decor') % Loading k space data for channel decorrelation
    data = complex(zeros(settings.parameters.pRO, settings.parameters.pPE, settings.parameters.pSlc, 'single'));
    for Slc = 1:settings.parameters.pSlc
        tmp = load([nameVol '/' nkspace '/Slice_' num2str(Slc) '/Channel_' num2str(Chunk) '.mat']);
        data(:,:,Slc) = tmp.data;
    end
end
