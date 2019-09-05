function settings=initializeProcess(settings, type)
% Function to initialize the parallel pool.
%
% ************************************************************************
% Version 1.0                                                  15.02.2019
% Falk Luesebrink              falk dot luesebrink at med dot ovgu dot de
% ************************************************************************

broken = false;
if strcmp(type, 'preproc')
    if settings.options.procParallel && ~settings.options.procMemory || settings.options.procParallel && settings.options.procMemory && settings.options.kspace.GRAPPA
        % Set maximum number of workers according to number of physical cores
        % (with an upper limit of 32 cores) and set maximum numbers of channels
        % to be processed in parallel (limited by available memory).
        if settings.options.autoSetup
            sliceMem = settings.parameters.nRO * settings.parameters.nPE * settings.parameters.nCh  *4 *2 / 10^9; % Size of slice times number of channels times 4 byte per element (single) times 2 (complex) divided by 10^9 to change unit to Gigabyte

            if settings.parameters.maxMemory == Inf
                if ispc
                    [~, sv] = memory;
                    settings.parameters.availableMem = sv.PhysicalMemory.Available / 10^9; % Change unit from Byte to GB
                else
                    [~, out]                         = system('vmstat -s -S M | grep "inactive memory"');
                    freeMem                          = sscanf(out,'%f free memory');
                    [~, out]                         = system('vmstat -s -S M | grep "free memory"');
                    inactiveMem                      = sscanf(out,'%f inactive memory');
                    settings.parameters.availableMem = (freeMem + inactiveMem) / 10^3; % Change unit from MB to GB
                end
            else
                settings.parameters.availableMem = settings.parameters.maxMemory;
            end

  if settings.options.memSafetyMargin
	  settings.parameters.availableMem = settings.parameters.availableMem * settings.parameters.safetyMargin;
  end

            if settings.parameters.preprocWorkers == Inf
                for i = 6:-1:0
                    settings.parameters.preprocWorkers = 2^i;
                    if settings.parameters.preprocWorkers <= feature('numcores') && settings.parameters.preprocWorkers <= settings.parameters.maxCPUs
                        for j = 1:32
                            settings.parameters.nChunks  = settings.parameters.preprocWorkers*j;
                            settings.parameters.stepSize = ceil(settings.parameters.pSlc / settings.parameters.nChunks);

                            if settings.options.kspace.GRAPPA
                                settings.parameters.estimatedSliceMem = 1.2 * sliceMem * (settings.parameters.nChunks * (settings.parameters.stepSize)) + 0.5;
                            else
                                settings.parameters.estimatedSliceMem = sliceMem * (settings.parameters.preprocWorkers * (settings.parameters.stepSize)) + 0.5;
                            end

                            if settings.parameters.availableMem > settings.parameters.estimatedSliceMem
                                broken = true;
                                break
                            end
                        end
                        if broken
                            break;
                        end
                    end
                end
            end
        else
            if isempty(settings.parameters.nChunks)
                settings.parameters.nChunks = settings.parameters.pSlc / settings.parameters.stepSize;
            elseif isempty(settings.parameters.stepSize)
                settings.parameters.stepSize = settings.parameters.pSlc / settings.parameters.nChunks;
            end
	    
	    sliceMem    = settings.parameters.nRO * settings.parameters.nPE * settings.parameters.nCh  *4 *2 / 1
	    0^9; % Size of slice times number of channels times 4 byte per element (single) times 2 (complex) divided by 10^9 to change unit to Gigabyte
	    if settings.options.kspace.GRAPPA
	    	settings.parameters.estimatedSliceMem = 1.2 * sliceMem * (settings.parameters.nChunks * (settings.parameters.stepSize)) + 0.5;
	    else
                settings.parameters.estimatedSliceMem = sliceMem * (settings.parameters.nChunks * (settings.parameters.stepSize)) + 0.5;
	    end
        end

        poolobj = gcp('nocreate');
        if isempty(poolobj)
            fprintf('| Starting parallel pool.\n');
            evalc('parpool(settings.parameters.preprocWorkers)');
            fprintf('| Connected to %i workers.\n', settings.parameters.preprocWorkers);
            fprintf('|\n');
        elseif poolobj.NumWorkers ~= settings.parameters.preprocWorkers
            fprintf('|\n');
            fprintf('| Parallel pool already running with different number of workers as requested.\n');
            fprintf('| Restarting parallel pool.\n');

            evalc('delete(gcp(''nocreate''))');
            evalc('parpool(settings.parameters.preprocWorkers)');

            fprintf('| Connected to %i workers.\n', settings.parameters.preprocWorkers);
            fprintf('|\n');
        end

        fprintf('| Total chunks to be processed: %i\n', settings.parameters.nChunks);
        fprintf('| Chunks processed in parallel: %i\n', settings.parameters.preprocWorkers);
        fprintf('| Each chunk consists of up to %i slices.\n', settings.parameters.stepSize);
	fprintf('| \n');
	fprintf('| Available memory: %.2f\n', settings.parameters.availableMem);
	if isfield(settings.parameters,'estimatedSliceMem');
		fprintf('| Estimated memory: %.2f\n', settings.parameters.estimatedSliceMem);
	end
	fprintf('| \n');
    else
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Setup chunks and stepSize
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if settings.options.autoSetup
            sliceMem    = settings.parameters.nRO * settings.parameters.nPE * settings.parameters.nCh  *4 *2 / 10^9; % Size of slice times number of channels times 4 byte per element (single) times 2 (complex) divided by 10^9 to change unit to Gigabyte
            
            if settings.parameters.maxMemory == Inf
                if ispc
                    [~, sv] = memory;
                    settings.parameters.availableMem = sv.PhysicalMemory.Available / 10^9; % Change unit from Byte to GB
                else
                    [~, out]                         = system('vmstat -s -S M | grep "inactive memory"');
                    freeMem                          = sscanf(out,'%f free memory');
                    [~, out]                         = system('vmstat -s -S M | grep "free memory"');
                    inactiveMem                      = sscanf(out,'%f inactive memory');
                    settings.parameters.availableMem = (freeMem + inactiveMem) / 10^3; % Change unit from MB to GB
                end
            else
                settings.parameters.availableMem = settings.parameters.maxMemory;
            end
            
            for j = 1:32
                settings.parameters.nChunks  = j;
                settings.parameters.stepSize = ceil(settings.parameters.pSlc / settings.parameters.nChunks);
                
                if settings.options.kspace.GRAPPA
                    settings.parameters.estimatedSliceMem = 1.2 * sliceMem * (settings.parameters.nChunks * (settings.parameters.stepSize)) + 0.5;
                else
                    settings.parameters.estimatedSliceMem = sliceMem * (settings.parameters.nChunks * (settings.parameters.stepSize)) + 0.5;
                end
                
                if settings.parameters.availableMem > settings.parameters.estimatedSliceMem
                    break
                end
            end
        else
            if isempty(settings.parameters.nChunks)
                settings.parameters.nChunks = settings.parameters.pSlc / settings.parameters.stepSize;
            elseif isempty(settings.parameters.stepSize)
                settings.parameters.stepSize = settings.parameters.pSlc / settings.parameters.nChunks;
            end
	end
    end
elseif strcmp(type, 'proc')
    if settings.options.procParallel
        channelMem  = settings.parameters.nRO * settings.parameters.nPE * settings.parameters.nSlc *4 *2 / 10^9; % Size of volume per channel times 4 byte per element (single) times 2 (complex) divided by 10^9 to change unit to Gigabyte
        % Set maximum number of workers according to number of physical cores
        % (with an upper limit of 32 cores) and set maximum numbers of channels
        % to be processed in parallel (limited by available memory).
        if settings.options.autoSetup
            if settings.parameters.maxMemory == Inf
                if ispc
                    [~, sv] = memory;
                    settings.parameters.availableMem = sv.PhysicalMemory.Available / 10^9; % Change unit from Byte to GB
                else
                    [~, out]                         = system('vmstat -s -S M | grep "inactive memory"');
                    freeMem                          = sscanf(out,'%f free memory');
                    [~, out]                         = system('vmstat -s -S M | grep "free memory"');
                    inactiveMem                      = sscanf(out,'%f inactive memory');
                    settings.parameters.availableMem = (freeMem + inactiveMem) / 10^3; % Change unit from MB to GB
                end
            else
                settings.parameters.availableMem = settings.parameters.maxMemory;
            end

            if settings.parameters.procWorkers == Inf
                for i = 6:-1:0
                    settings.parameters.procWorkers = 2^i;
                    if settings.parameters.procWorkers <= feature('numcores') && settings.parameters.procWorkers <= settings.parameters.maxCPUs 
                        settings.parameters.procWorkers = settings.parameters.procWorkers;
                        if settings.options.kspace.denoise
                            settings.parameters.estimatedChannelMem = (channelMem * settings.parameters.procWorkers * 6);
			else
                            settings.parameters.estimatedChannelMem = (channelMem * settings.parameters.procWorkers * 2);
                        end
                        if settings.parameters.availableMem > settings.parameters.estimatedChannelMem
                            break;
                        end
                    end
                end
            end
	else
		if settings.options.kspace.denoise
			settings.parameters.estimatedChannelMem = channelMem * settings.parameters.procWorkers * 2;
		else
			settings.parameters.estimatedChannelMem = channelMem * settings.parameters.procWorkers * 6;
		end
        end

        if settings.parameters.procWorkers > settings.parameters.nCh
            settings.parameters.procWorkers = settings.parameters.nCh;
        end

        poolobj = gcp('nocreate');
        if isempty(poolobj)
            fprintf('| Starting parallel pool.\n');
            evalc('parpool(settings.parameters.procWorkers)');
            fprintf('| Connected to %i workers.\n', settings.parameters.procWorkers);
            fprintf('|\n');
        elseif poolobj.NumWorkers ~= settings.parameters.procWorkers
            fprintf('| \n');
            fprintf('| Parallel pool already running with different number of workers as requested.\n');
            fprintf('| Restarting parallel pool.\n');

            evalc('delete(gcp(''nocreate''))');
            evalc('parpool(settings.parameters.procWorkers)');

            fprintf('| Connected to %i workers.\n', settings.parameters.procWorkers);
            fprintf('|\n');
        end

        fprintf('| Total chunks to be processed: %i\n', settings.parameters.nCh);
        fprintf('| Chunks processed in parallel: %i\n', settings.parameters.procWorkers);
        fprintf('| Each chunk contains data of one channel.\n');
	fprintf('| \n');
	fprintf('| Available memory: %.2f\n', settings.parameters.availableMem);
	%fprintf('| Estimated memory: %.2f\n', settings.parameters.estimatedChannelMem);
        fprintf('|\n');
    end
elseif strcmp(type, 'combine')
    if settings.options.procParallel
        % Set maximum number of workers according to number of physical cores
        % (with an upper limit of 32 cores) and set maximum numbers of channels
        % to be processed in parallel (limited by available memory).
        if settings.options.autoSetup
            sliceMem    = settings.parameters.nRO * settings.parameters.nPE * settings.parameters.nCh  *4 *2 / 10^9; % Size of slice times number of channels times 4 byte per element (single) times 2 (complex) divided by 10^9 to change unit to Gigabyte

            if settings.parameters.maxMemory == Inf
                if ispc
                    [~, sv] = memory;
                    settings.parameters.availableMem = sv.PhysicalMemory.Available / 10^9; % Change unit from Byte to GB
                else
                    [~, out]                         = system('vmstat -s -S M | grep "inactive memory"');
                    freeMem                          = sscanf(out,'%f free memory');
                    [~, out]                         = system('vmstat -s -S M | grep "free memory"');
                    inactiveMem                      = sscanf(out,'%f inactive memory');
                    settings.parameters.availableMem = (freeMem + inactiveMem) / 10^3; % Change unit from MB to GB
                end
            else
                settings.parameters.availableMem = settings.parameters.maxMemory;
            end

            if settings.parameters.combineWorkers == Inf
                for j = 1:32
                    for i = 6:-1:0
                        settings.parameters.combineWorkers = 2^i;
                        if settings.parameters.combineWorkers <= feature('numcores') && settings.parameters.combineWorkers <= settings.parameters.maxCPUs
                            settings.parameters.nChunks  = settings.parameters.combineWorkers*j;
                            settings.parameters.stepSize = ceil(settings.parameters.nSlc / settings.parameters.nChunks);
                            settings.parameters.estimatedSliceMem = 2 * sliceMem * (settings.parameters.combineWorkers * (settings.parameters.stepSize));
                            if settings.parameters.availableMem > settings.parameters.estimatedSliceMem
                                broken = true;
                                break
                            end
                        end
                    end
                    if broken
                        break;
                    end
                end
            end
        else
            if isempty(settings.parameters.nChunks)
                settings.parameters.nChunks = settings.parameters.nSlc / settings.parameters.stepSize;
            elseif isempty(settings.parameters.stepSize)
                settings.parameters.stepSize = settings.parameters.nSlc / settings.parameters.nChunks;
            end
	    sliceMem    = settings.parameters.nRO * settings.parameters.nPE * settings.parameters.nCh  *4 *2 / 10^9; % Size of slice times number of channels times 4 byte per element (single) times 2 (complex) divided by 10^9 to change unit to Gigabyte
	    settings.parameters.estimatedSliceMem = 2 * sliceMem * (settings.parameters.combineWorkers * (settings.parameters.stepSize));
        end

        poolobj = gcp('nocreate');
        if isempty(poolobj)
            fprintf('| Starting parallel pool.\n');
            evalc('parpool(settings.parameters.combineWorkers)');
            fprintf('| Connected to %i workers.\n', settings.parameters.combineWorkers);
            fprintf('|\n');
        elseif poolobj.NumWorkers ~= settings.parameters.combineWorkers
            fprintf('|\n');
            fprintf('| Parallel pool already running with different number of workers as requested.\n');
            fprintf('| Restarting parallel pool.\n');

            evalc('delete(gcp(''nocreate''))');
            evalc('parpool(settings.parameters.combineWorkers)');

            fprintf('| Connected to %i workers.\n', settings.parameters.combineWorkers);
            fprintf('|\n');
        end

        fprintf('| Total chunks to be processed: %i\n', settings.parameters.nChunks);
        fprintf('| Chunks processed in parallel: %i\n', settings.parameters.combineWorkers);
        fprintf('| Each chunk consists of up to %i slices.\n', settings.parameters.stepSize);
        fprintf('| \n');
	fprintf('| Available memory: %.2f\n', settings.parameters.availableMem);
	%fprintf('| Estimated memory: %.2f\n', settings.parameters.estimatedSliceMem);
	fprintf('| \n');
    else
        if settings.options.autoSetup
            sliceMem    = settings.parameters.nRO * settings.parameters.nPE * settings.parameters.nCh  *4 *2 / 10^9; % Size of slice times number of channels times 4 byte per element (single) times 2 (complex) divided by 10^9 to change unit to Gigabyte

            if settings.parameters.maxMemory == Inf
                if ispc
                    [~, sv] = memory;
                    settings.parameters.availableMem = sv.PhysicalMemory.Available / 10^9; % Change unit from Byte to GB
                else
                    [~, out]                         = system('vmstat -s -S M | grep "inactive memory"');
                    freeMem                          = sscanf(out,'%f free memory');
                    [~, out]                         = system('vmstat -s -S M | grep "free memory"');
                    inactiveMem                      = sscanf(out,'%f inactive memory');
                    settings.parameters.availableMem = (freeMem + inactiveMem) / 10^3; % Change unit from MB to GB
                end
            else
                settings.parameters.availableMem = settings.parameters.maxMemory;
            end

            for j = 1:32
                settings.parameters.nChunks  = j;
                settings.parameters.stepSize = ceil(settings.parameters.nSlc / settings.parameters.nChunks);
                settings.parameters.estimatedSliceMem = 2 * sliceMem * (settings.parameters.combineWorkers * (settings.parameters.stepSize));

                if settings.parameters.availableMem > settings.parameters.estimatedSliceMem
                    break
                end
            end
        else
            if isempty(settings.parameters.nChunks)
                settings.parameters.nChunks = settings.parameters.nSlc / settings.parameters.stepSize;
            elseif isempty(settings.parameters.stepSize)
                settings.parameters.stepSize = settings.parameters.nSlc / settings.parameters.nChunks;
            end
        end
    end
end
	
if settings.parameters.preprocWorkers == 1 && settings.parameters.procWorkers == 1 && settings.options.procParallel
	settings.options.procParallel = false;
end
