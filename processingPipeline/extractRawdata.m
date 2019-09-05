function kspace=extractRawdata(name, settings, twix_obj, Chunk)
% Extracting raw data, apply decorrelation (pre-whitening) and write it to
% disk. 
%
% ************************************************************************
% Version 1.0                                                  15.02.2019
% Falk Luesebrink              falk dot luesebrink at med dot ovgu dot de
% ************************************************************************
if ~settings.options.procParallel || settings.options.procMemory && settings.options.procParallel && ~settings.options.kspace.GRAPPA
	kspace = twix_obj.image{''};
	kspace = permute(kspace, [2 1 3 4]);
    Index  = 1:settings.parameters.pSlc;
    Chunk  = 1;
else
	Index = ((Chunk-1)*settings.parameters.stepSize)+1 : Chunk*settings.parameters.stepSize;

	% If index requests a slice beyond the maximum number of slices, remove that element
	if Index(end) > settings.parameters.pSlc
		stepSize = settings.parameters.stepSize - (Index(end) - settings.parameters.pSlc); 
		Index = Index(1:stepSize);
    end
    
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% Extract k space data from raw data file.
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	kspace = complex(zeros(settings.parameters.pRO, settings.parameters.nCh, settings.parameters.pPE, settings.parameters.stepSize, 'single'));
	
	%MQ
    % store data of first average
    if twix_obj.image.NAve > 1 && settings.options.kspace.variableAveraging 
        data = squeeze(twix_obj.image(:,:,:,:,1,1)); 
    end
	
	for Line = 1:settings.parameters.pPE
		%MQ: average the data, if variable averaging is enabled and the
        if twix_obj.image.NAve > 1 && settings.options.kspace.variableAveraging 
            
            kspace(:,:,Line,Index-(Chunk-1)*settings.parameters.stepSize) = data(:,:,Line,Index);
            
            % compute the AverageFactor
            AverageFactors = getAverageFactors(twix_obj, Line, (Index(1):Index(end)));

            % add the data from the remaining averages
            %for Average = 2:twix_obj.image.NAve
            if max(AverageFactors)>1
                % loop to 
                for Average = 2:max(AverageFactors) 
                    % load the additional data
                    tmp_ave = squeeze(twix_obj.image(:,:,:,Index,1,Average)); 
                    % add the average data to the data
                    kspace(:,:,Line,Index-(Chunk-1)*settings.parameters.stepSize) = tmp_ave(:,:,Line,:) + ...
                        kspace(:,:,Line,Index-(Chunk-1)*settings.parameters.stepSize); 
                end

                % correct for the AverageFactor
                for Partition=1:length(Index)
                    kspace(:,:,Line,Partition) = kspace(:,:,Line,Partition) * (1/squeeze(AverageFactors(1,Partition)));
                end
            end
        else
            data = twix_obj.image(:,:,Line,:);
            kspace(:,:,Line,Index-(Chunk-1)*settings.parameters.stepSize) = data(:,:,1,Index);
        end

    % only if acquisition weighting is enabled
        if settings.options.kspace.acquisitionWeighting
            % set the amplitude of the acquisition weighting
            if twix_obj.image.NAve > 1
                AW_amplitude = twix_obj.image.NAve;
            else
                AW_amplitude = 1;
            end
            
            AcqWeightingFactors = getAcqWeightingFactors(settings, Line, (Index(1):Index(end)), ...
                AW_amplitude);

            for Partition=1:length(Index)
                kspace(:,:,Line,Partition) = kspace(:,:,Line,Partition) * AcqWeightingFactors(1, Partition);
            end
        end
        %MQ end
	end
	kspace = permute(kspace, [2 1 3 4]);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Decorrelate channels (optional)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if settings.options.kspace.decorrelateChannels && ~isscalar(settings.parameters.W)
    kspace = reshape(settings.parameters.W * kspace(:,:), size(kspace));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Apply GRAPPA reconstruction (if needed)
%
% Does not seem to work properly. And ghosting is still
% present. I haven't found a solution and probably can
% fixed by using different code.
%
% Should be moved to 'preprocRawdata.m' to allow for GRAPPA
% reconstruction in other dimensions. This makes it
% necessary to move the save function as well.
%
% Cannot be used with denoising in its current state as
% noise distribution will be changed due to the interpolation
% with the GRAPPA kernel.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if settings.options.kspace.GRAPPA
%     tmp = complex(zeros(32, 224, 224, 160, 'single'));
%     tmp(:,:,2:223,:) = kspace;
%     kspace = tmp;
%     clear tmp

    kspace = permute(kspace, [2 3 1 4]);
	refscan = twix_obj.refscan{''};
	refscan = permute(refscan, [1 3 2 4]);
    refscan = refscan(:,:,:,Index);
	kSize = [3,3];
	for Slice = 1:length(Index)
		kspace(:,:,:,Slice) = applyGRAPPA(kspace(:,:,:,Slice), refscan(:,:,:,Slice), kSize, 0.1);
	end
	kspace = permute(kspace, [3 1 2 4]);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Save to disk
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if settings.options.saveIntermediate && settings.options.procMemory
	saveData(name, settings, kspace, Index, Chunk, 'preproc')
% Is this acutually needed? And probably saveIntermediate can be used if
% procMemory is used anyways, making the && in the statement above obsolet.
% elseif settings.options.saveIntermediate 
%     saveData(name, settings, kspace, [], [], 'kspace_decor')
elseif ~settings.options.procMemory
	saveData(name, settings, kspace, Index, Chunk, 'preproc')
	kspace = 0;
end