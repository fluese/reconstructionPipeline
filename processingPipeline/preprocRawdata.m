function [data,settings]=preprocRawdata(name, settings, twix_obj)
% Function to reconstruct raw MRI data into x space.
%
% ************************************************************************
% Version 1.0                                                  15.02.2019
% Falk Luesebrink              falk dot luesebrink at med dot ovgu dot de
% ************************************************************************

fprintf('===============================================================\n');
fprintf('| Pre-processing of raw data started at %s\n', datestr(datetime('now')));
fprintf('===============================================================\n');
tstart=tic;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Check if preprocessed data already exists
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
settings = statusData(name, settings, 'kspace');

if ~settings.parameters.kspace.processed
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Display what processes are run here!!!!!!!!!!!!!!!!!!!!!
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    fprintf('| Preprocessing k space data...\n');
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Initialize Process
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    settings=initializeProcess(settings, 'preproc');

    if settings.options.procMemory
        if settings.options.kspace.GRAPPA && settings.options.procParallel
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Parallel stream (works for GRAPPA, but not in general)
        %
        % It seems to be faster to run this not in parallel. Probably
        % due to the cell structure. If the data is processed in memory
        % it doesn't make much sense to run it in parallel in its
        % current state.
        %
        % Find a better solution than using the cells and it will
        % run faster!
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            tmp=cell(1,settings.parameters.nChunks);
            fprintf(['| ' repmat('.',1,settings.parameters.nChunks) '\n| \n']);
            parfor Chunk = 1:settings.parameters.nChunks
                fprintf('\b|\n');
              
                tmp{1,Chunk} = extractRawdata(name, settings, twix_obj, Chunk);
                tmp{1,Chunk} = permute(tmp{1,Chunk}, [1 4 2 3]);
            end
            data = cell2mat(tmp);
            data = permute(data,[1 3 4 2]);
            data = data(:,:,:,1:settings.parameters.pSlc);
        else
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Single CPU stream
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            disp(['| Processing slices 1 to ' num2str(settings.parameters.pSlc)]);
            data = extractRawdata(name, settings, twix_obj);
        end
    else
        if settings.options.procParallel
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Parallel stream
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            fprintf(['| ' repmat('.',1,settings.parameters.nChunks) '\n| \n']);
            parfor Chunk=1:settings.parameters.nChunks
                fprintf('\b|\n');
                extractRawdata(name, settings, twix_obj, Chunk);
            end
        else
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Single CPU stream
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            for Chunk = 1 : settings.parameters.nChunks
                if Chunk*settings.parameters.stepSize >= settings.parameters.pSlc
                    disp(['| Processing slices ' num2str(((Chunk-1)*settings.parameters.stepSize)+1) ' to ' num2str(settings.parameters.pSlc)]);
                else
                    disp(['| Processing slices ' num2str(((Chunk-1)*settings.parameters.stepSize)+1) ' to ' num2str(Chunk*settings.parameters.stepSize)]);
                end

                extractRawdata(name, settings, twix_obj, Chunk);
            end
        end
        data = 0;
    end
else
    fprintf('| Raw data already extracted.\n');
    fprintf('===============================================================\n');
    data = 0;
    return;
end
tend=toc(tstart);
fprintf('| \n');
disp(['| Done pre-processing k space data. Runtime: ' datestr(tend/(24*60*60), 'DD:HH:MM:SS.FFF')]);
fprintf('===============================================================\n');    