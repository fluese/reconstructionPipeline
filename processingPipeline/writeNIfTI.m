function writeNIfTI(name, settings, imspace)
% Function to write NIfTI files.
%
% ************************************************************************
% Version 1.0                                                  15.02.2019
% Falk Luesebrink              falk dot luesebrink at med dot ovgu dot de
% ************************************************************************

fprintf('===============================================================\n');
fprintf('| Generating NIfTI started at %s \n', datestr(datetime('now')));
fprintf('===============================================================\n');

settings=statusData(name, settings, 'NIfTI');

if settings.options.output.writeNIfTI && ~settings.parameters.NIfTI.processed
    if ~settings.options.procMemory || ~(size(imspace,1) > 1)
        % Setup variable
        imspace = complex(zeros(settings.parameters.nRO, settings.parameters.nPE, settings.parameters.nSlc, 'single'));

        fprintf('| Loading reconstructed data...\n');
        for k=1:settings.parameters.nSlc % use parload
            tmp=load([name.namePath '/' name.imspace '/Slice_' num2str(k) '.mat']);
            imspace(:,:,k)=tmp.data;
        end
    end
    
    % Preparing orientation of data.
    imspace = permute(imspace, [2 1 3]);
    imspace = imspace(end:-1:1,end:-1:1,:);
    
    % Get magnitude image and scale it.
    vol = abs(imspace);
    vol = ( (vol - min(vol(:)) ) / (max(vol(:)) - min(vol(:)) ) ) .* settings.parameters.scalefac;
    
    % Create NIfTI file
    data=make_nii(vol, settings.parameters.Resolution, [], 16);
    
    % Set NIfTI magic number
    data.hdr.hist.magic = 'n+1';
    
    % Set unit of voxel spacing to mm
    data.hdr.dime.xyzt_units = 2;
    
    % Change type of orientation method (qform: quanternion based; sform:
    % affine transformation). Compare nifti header information for information
    % on different methods.
    data.hdr.hist.qform_code = 1; % Method 2: qform
    % data.hdr.hist.sform_code = 1; % Method 3: sform
    
    % Quaternion (Method 2)
    data.hdr.hist.quatern_b =    settings.parameters.Quaternion(1,3);
    data.hdr.hist.quatern_c = -1*settings.parameters.Quaternion(1,2);
    data.hdr.hist.quatern_d = -1*settings.parameters.Quaternion(1,1);
%     data.hdr.hist.quatern_b = -1 * settings.parameters.Quaternion(1,3);
%     data.hdr.hist.quatern_c =  1 * settings.parameters.Quaternion(1,2);
%     data.hdr.hist.quatern_d = -1 * settings.parameters.Quaternion(1,1);
    data.hdr.hist.qoffset_x = 0; %settings.parameters.FoV(3)/2 + -1 * (settings.parameters.Pos(1) + settings.parameters.Resolution(1)/2);
    data.hdr.hist.qoffset_y = 0; %settings.parameters.FoV(2)/2 + -1 * settings.parameters.Pos(2);
    data.hdr.hist.qoffset_z = 0; %settings.parameters.FoV(1)/2 +      settings.parameters.Pos(3);
    data.hdr.dime.pixdim(1) = -1;
    
    % % Affine transformation part (Method 3)
    data.hdr.hist.srow_x =  settings.parameters.R(1,1:4);
    data.hdr.hist.srow_y =  settings.parameters.R(2,1:4);
    data.hdr.hist.srow_z =  settings.parameters.R(3,1:4);
    
    % Save NIfTI to disk
    if ~settings.options.output.writePhase
        fprintf('| Writing NIfTI file...\n');
        save_nii(data, [name.namePath '/' name.nameFile '.nii'])
        if settings.options.output.gzip
            gzip([name.namePath '/' name.nameFile '.nii']);
            delete([name.namePath '/' name.nameFile '.nii']);
        end
    elseif settings.options.output.writePhase && strcmp(settings.options.xspace.combineChannels,'adaptComb')
        fprintf('| Writing magnitude NIfTI file...\n');
        save_nii(data, [name.namePath '/' name.nameFile '.nii'])
        if settings.options.output.gzip
            gzip([name.namePath '/' name.nameFile '.nii']);
            delete([name.namePath '/' name.nameFile '.nii']);
        end
        
        fprintf('| Writing phase NIfTI file...\n');
        if settings.options.output.gzip
            tmp = load_untouch_nii([name.namePath '/' name.nameFile '.nii.gz']);
        else
            tmp = load_untouch_nii([name.namePath '/' name.nameFile '.nii']);
        end
        tmp = rmfield(tmp,'untouch');
        tmp.img = angle(imspace);
        
        save_nii(tmp, [name.namePath '/' name.nameFile '_phase.nii']);
        
        if settings.options.output.gzip
            data=[name.namePath '/' name.nameFile '_phase.nii'];
            gzip(data)
            delete(data)
        end
    elseif settings.options.output.writePhase && ~strcmp(settings.options.xspace.combineChannels,'adaptComb')
        fprintf('| Writing magnitude NIfTI file only.\n');
        fprintf('| Phase information is omitted in channel combination by sum of squares.\n');
        save_nii(data, [name.namePath '/' name.nameFile '.nii'])
        if settings.options.output.gzip
            gzip([name.namePath '/' name.nameFile '.nii']);
            delete([name.namePath '/' name.nameFile '.nii']);
        end
    end
    
    if settings.options.output.writeComplex
        % As complex data cannot really be handled by neuroimaging tools we
        % transform the 3D image into a 4D time series. The real part will
        % be the first time point, the imaginary the second one.
        % This way transforms (e.g. generated by ANTs) can be applied to
        % the dataset. Afterwards the split complex number can be merged
        % again.
        fprintf('| Writing complex valued NIfTI file additionally...\n');
        if settings.options.output.gzip
            tmp = load_untouch_nii([name.namePath '/' name.nameFile '.nii.gz']);
        else
            tmp = load_untouch_nii([name.namePath '/' name.nameFile '.nii']);
        end
        tmp = rmfield(tmp,'untouch');
        tmp.img = [];
        
        tmp.hdr.dime.datatype = 16;
        tmp.hdr.dime.bitpix = 32;
        
        tmp.hdr.dime.dim(1) = 4;
        tmp.hdr.dime.dim(5) = 2;
        tmp.img(:,:,:,1) = real(imspace);
        tmp.img(:,:,:,2) = imag(imspace);
        save_nii(tmp, [name.namePath '/' name.nameFile '_complex.nii']);
        
        if settings.options.output.gzip
            data=[name.namePath '/' name.nameFile '_complex.nii'];
            gzip(data)
            delete(data)
        end
        
        % Write true complex NIfTI
        if settings.options.output.trueComplex
            tmp.hdr.dime.datatype = 32;
            tmp.hdr.dime.bitpix = 64;
            
            tmp.img = imspace;
            save_nii(tmp, [name.namePath '/' name.nameFile '_trueComplex.nii']);
            if settings.options.output.gzip
                data=[name.namePath '/' name.nameFile '_trueComplex.nii'];
                gzip(data)
                delete(data)
            end
        end
    end
    
    if settings.options.output.refData
        if ~exist(settings.parameters.output.ref_nifti, 'file')
            fprintf('| The specified reference image does not exist. Skipping to write data based on reference.\n');
        else
            fprintf('| Reading header information of reference NIfTI to get orientation information.\n');
            [~, ~, ext] = fileparts(settings.parameters.output.ref_nifti);
            if strcmp(ext, '.gz')
                gunzip(settings.parameters.output.ref_nifti);
                HeaderInfo = spm_vol(settings.parameters.output.ref_nifti(1:end-3));
            else
                HeaderInfo = spm_vol(settings.parameters.output.ref_nifti);
            end
            
            ref = spm_read_vols(HeaderInfo);
	    scale = max(ref(:)); clear ref;
            
            vol=circshift(vol,[1,0,-1]); % It seems as if there is a shift between data reconstructed at the scanner and my reconstruction pipeline. Reason?

	    vol = vol(:,end:-1:1,end:-1:1);
	    vol = vol ./ max(vol(:)) .* scale;
	                
            fprintf('| Writing NIfTI file based on reference NIfTI...\n');
            HeaderInfo.fname=([name.namePath '/' name.nameFile '_withRef.nii']);
            HeaderInfo.private.dat.fname=HeaderInfo.fname;
            HeaderInfo.dt=[512 0];
            spm_write_vol(HeaderInfo,vol);
            
            if settings.options.output.gzip
                gzip(HeaderInfo.fname)
                delete(HeaderInfo.fname)
            end
        end
    end
    
    fprintf('| \n');
    fprintf('| Done writing NIfTI files to disk.\n')
    fprintf('===============================================================\n');
elseif settings.parameters.NIfTI.processed
    fprintf('| NIfTI file already exists.\n')
    fprintf('===============================================================\n');
end
