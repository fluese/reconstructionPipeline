function applyDistortionCorrection(name, settings)
% Function to apply 3D distortion correction retrospectively using code of
% Uten Yarach.
%
% ************************************************************************
% Version 1.0                                                  15.02.2019
% Falk Luesebrink              falk dot luesebrink at med dot ovgu dot de
% ************************************************************************

if settings.options.post.distortionCorrection
    if settings.options.output.refData
        fprintf('===============================================================\n');
        fprintf('| Conducting distortion correction. Started at %s\n', datestr(datetime('now')));
        fprintf('===============================================================\n');
        
        if settings.options.output.gzip
            refName = [name.namePath '/' name.nameFile '_withRef.nii.gz'];
        else
            refName = [name.namePath '/' name.nameFile '_withRef.nii'];
        end
        
        % Read data and orientation
        [img, RAS]= mris_read_nii(refName);
        
        % apply DiCo to Image
        dicoimg = CorrectDistortion(img, RAS);
        
        % Write NIfTI file
        fprintf('| Writing NIfTI file...\n')
        [~,~,ext]=fileparts(refName);
        if strcmp(ext, '.gz')
            tmp=refName(1:end-3);
            HeaderInfo=spm_vol_nifti(tmp);
            HeaderInfo.fname = [tmp(1:end-4) '_DiCo.nii'];
        else
            HeaderInfo=spm_vol_nifti(refName);
            HeaderInfo.fname = [refName(1:end-4) '_DiCo.nii'];
        end
        HeaderInfo.private.dat.fname = HeaderInfo.fname;
        
        spm_write_vol(HeaderInfo,dicoimg);
        if settings.options.output.gzip
            gzip(HeaderInfo.fname);
            delete(HeaderInfo.fname);
        end
        fprintf('| Done.\n')
        fprintf('===============================================================\n\n');
    else
        fprintf('===============================================================\n');
        fprintf('| Conducting distortion correction. Started at %s\n', datestr(datetime('now')));
        fprintf('===============================================================\n');
        fprintf('| Currently the origin is set to the center of the image. Please use a reference image to conduct distortion correction.\n');
        fprintf('===============================================================\n\n');
    end
end