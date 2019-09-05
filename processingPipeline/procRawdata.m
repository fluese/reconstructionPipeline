function [data, settings]=procRawdata(name, settings, data)
% Function to reconstruct raw MRI data into x space.
%
% ************************************************************************
% Version 1.0                                                  15.02.2019
% Falk Luesebrink              falk dot luesebrink at med dot ovgu dot de
% ************************************************************************

fprintf('\n')
fprintf('===============================================================\n');
fprintf('| Processing of raw data started at %s\n', datestr(datetime('now')));
fprintf('===============================================================\n');
tstart=tic;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Check if processed data already exists
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
settings = statusData(name, settings, 'xspace');

if ~settings.parameters.xspace.processed
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Display what processes are run here!!!!!!!!!!!!!!!!!!!!!
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    fprintf('| Processing k space data...\n');
            
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Initialize Process
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    settings=initializeProcess(settings, 'proc');
    
    if settings.options.procParallel
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Parallel stream
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            fprintf(['| ' repmat('.',1,settings.parameters.nCh) '\n| \n']);
            if settings.options.procMemory
                if data == 0
                    data = complex(zeros(settings.parameters.nCh, settings.parameters.pRO, settings.parameters.pPE, settings.parameters.pSlc, 'single'));
                    for Channel = 1:settings.parameters.nCh
                        data(Channel,:,:,:) = loadData(name, settings, [], Channel, 'kspace_decor');
                    end
                end
                
                if settings.options.kspace.partialFourier
                    % There has to be a better way to do this!
                    tmp = complex(zeros(settings.parameters.nCh, settings.parameters.nRO, settings.parameters.nPE, settings.parameters.nSlc, 'single'));
                    parfor Channel = 1 : settings.parameters.nCh
                        fprintf('\b|\n');
                        tmp(Channel,:,:,:)=reconstructData(name, settings, Channel, squeeze(data(Channel,:,:,:)));
                    end
                    data = tmp; clear tmp
                else
                    parfor Channel = 1 : settings.parameters.nCh
                        fprintf('\b|\n');
                        data(Channel,:,:,:)=reconstructData(name, settings, Channel, squeeze(data(Channel,:,:,:)));
                    end
                end
            else
                parfor Channel = 1 : settings.parameters.nCh
                    fprintf('\b|\n');
                    reconstructData(name, settings, Channel);
                end
            end

            if ~settings.options.procMemory
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                % Reformat data into single files per slice
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                data = 0;
                fprintf('| \n');
                fprintf('| Separated channels will be written in a single file per slice.\n');
                saveData(name, settings, [], [], [], 'sepChannels');
            elseif settings.options.saveIntermediate
                saveData(name, settings, data, [], [], 'saveIntermediate');
            end
    else
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Single CPU stream
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if settings.options.procMemory && settings.options.saveIntermediate
            data=complex(zeros(settings.parameters.nCh, settings.parameters.pRO, settings.parameters.pPE, settings.parameters.pSlc, 'single'));
            for Channel = 1:settings.parameters.nCh
                data(Channel,:,:,:) = loadData(name, settings, [], Channel, 'kspace_decor');
            end
        end
        
        if settings.options.procMemory
            % Is there a better to do this without a tmp variable?
            tmp = complex(zeros(settings.parameters.nCh, settings.parameters.nRO, settings.parameters.nPE, settings.parameters.nSlc, 'single'));
            disp(['| Processing channel 1 to ' num2str(settings.parameters.nCh)])
            for Channel = 1:settings.parameters.nCh
                tmp(Channel,:,:,:)=reconstructData(name, settings, Channel, squeeze(data(Channel,:,:,:)));
            end
            data = tmp; clear tmp
        else
            for Channel = 1 : settings.parameters.nCh
                disp(['| Processing channel ' num2str(Channel)])
                reconstructData(name, settings, Channel);
            end
        end
        
        if ~settings.options.procMemory
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Reformat data into single files per slice
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            data = 0;
            fprintf('| \n');
            fprintf('| Separated channels will be written in a single file per slice.\n');
            saveData(name, settings, [], [], [], 'sepChannels');
        elseif settings.options.saveIntermediate
            saveData(name, settings, data, [], [], 'saveIntermediate');
        end
    end
else
    fprintf('| Data already reconstructed into x space.\n');
    fprintf('===============================================================\n');
    data = 0;
    return;
end
tend=toc(tstart);
fprintf('| \n');
disp(['| Done reconstructing k space into x space. Runtime: ' datestr(tend/(24*60*60), 'DD:HH:MM:SS.FFF')]);
fprintf('===============================================================\n');