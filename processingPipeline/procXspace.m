function [imspace, settings]=procXspace(name, settings, xspace)
% Function to combine MRI data from x space into image space.
%
% ************************************************************************
% Version 1.0                                                  15.02.2019
% Falk Luesebrink              falk dot luesebrink at med dot ovgu dot de
% ************************************************************************

fprintf('\n')
fprintf('===============================================================\n');
fprintf('| Processing of data in x space started at %s\n', datestr(datetime('now')));
fprintf('===============================================================\n');
tstart=tic;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Check if data already exists
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
settings=statusData(name, settings, 'imspace');

if ~settings.parameters.imspace.processed
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Display what processes are run here!!!!!!!!!!!!!!!!!!!!!!
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    fprintf('| Combing channels.\n');
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Initialize Process
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    settings=initializeProcess(settings, 'combine');
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Retrieve channel with highest intensity for phase correction
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    settings=getMaxChannel(name, settings, xspace);
    
    if settings.options.procParallel
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Parallel stream
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        fprintf(['| ' repmat('.',1,settings.parameters.nChunks) '\n| \n']);
        if settings.options.procMemory
            if settings.options.procMemory && settings.options.saveIntermediate
                Index = 1 : settings.parameters.nSlc;
                xspace = loadData(name, settings, Index, [], 'saveIntermediate');
            end
            
            tmp=cell(1,settings.parameters.nChunks);
            for Chunk = 1:settings.parameters.nChunks
                Index = ((Chunk-1)*settings.parameters.stepSize)+1 : Chunk*settings.parameters.stepSize;
                
                if Index(end) > settings.parameters.nSlc
                    entriesToRemove = Index > settings.parameters.nSlc;
                    Index(entriesToRemove) = [];
                end
                
                tmp{1,Chunk}=xspace(:,:,:,Index);
            end
            
            parfor Chunk=1:settings.parameters.nChunks
                fprintf('\b|\n');
                tmp{1,Chunk} = combineChannels(name, settings, Chunk, tmp{1,Chunk});
                tmp{1,Chunk} = permute(tmp{1,Chunk}, [1 3 2]);
            end
            
            imspace = cell2mat(tmp);
            imspace = permute(imspace,[1 3 2]);
        else
            parfor Chunk=1:settings.parameters.nChunks
                fprintf('\b|\n');
                combineChannels(name, settings, Chunk);
            end
            imspace = 0;
        end
    else
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Single CPU stream
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if settings.options.procMemory
            disp(['| Processing slices 1 to ' num2str(settings.parameters.nSlc)]);
            if xspace == 0
                Index = 1 : settings.parameters.nSlc;
                xspace = loadData(name, settings, Index, [], 'saveIntermediate');
            end
            settings.parameters.nChunks  = 1;
            settings.parameters.stepSize = settings.parameters.nSlc;
            imspace=combineChannels(name, settings, 1, xspace);
        else
            for Chunk=1:settings.parameters.nChunks
                if ((Chunk-1)*settings.parameters.stepSize)+1 > settings.parameters.nSlc
                    break;
                elseif Chunk*settings.parameters.stepSize > settings.parameters.nSlc
                    disp(['| Processing slices ' num2str(((Chunk-1)*settings.parameters.stepSize)+1) ' to ' num2str(settings.parameters.nSlc)]);
                    combineChannels(name, settings, Chunk);
                else
                    disp(['| Processing slices ' num2str(((Chunk-1)*settings.parameters.stepSize)+1) ' to ' num2str(Chunk*settings.parameters.stepSize)]);
                    combineChannels(name, settings, Chunk);
                end
            end
            imspace = 0;
        end
    end
else
    fprintf('| Data already in image space.\n');
    fprintf('===============================================================\n\n');
    imspace = 0;
    return;
end
tend=toc(tstart);
fprintf('| \n');
disp(['| Done reconstruction x space into image space. Runtime: ' datestr(tend/(24*60*60), 'DD:HH:MM:SS.FFF')]);
fprintf('===============================================================\n\n');