function applyBiasCorrection(name, settings)
% Apply bias field correction
%
% ************************************************************************
% Version 1.0                                                  15.02.2019
% Falk Luesebrink              falk dot luesebrink at med dot ovgu dot de
% ************************************************************************


if settings.options.post.biasfield
    if ~settings.parameters.SPM
        fprintf('\n===============================================================\n');
        fprintf('| You need to download SPM12 and add it to your Matlab path in \n');
        fprintf('| order to run the bias field correction. \n');
        fprintf('===============================================================\n');
    else
        if settings.options.post.segmentation
            fprintf('\n===============================================================\n');
            fprintf('| Conducting bias field correction and segmentation. Started at %s\n', datestr(datetime('now')));
            fprintf('===============================================================\n');
        else
            fprintf('\n===============================================================\n');
            fprintf('| Conducting bias field correction. Started at %s\n', datestr(datetime('now')));
            fprintf('===============================================================\n');
        end

        [fpath,fname,ext]=fileparts(name.NIfTI);
        if strcmp(ext, '.gz')
            gunzip(name.NIfTI);
            [fpath, fname, ext]=fileparts(name.NIfTI(1:end-3));
        end
	biasCorrection(settings, [fpath '/' fname ext]);

        if settings.options.post.mask
            % Creating brain mask by combining the segmentations, thresholding
            % it at greater than 0.5, and applying morphological operators to
            % refine it.
            fprintf('|\n')
            fprintf('| Creating brain mask from segmentation.\n')
            fprintf('|\n')

            HeaderInfo = spm_vol_nifti([fpath '/c1' fname ext]);
            i1 = spm_read_vols(HeaderInfo);
            HeaderInfo = spm_vol_nifti([fpath '/c2' fname ext]);
            i2 = spm_read_vols(HeaderInfo);
            HeaderInfo = spm_vol_nifti([fpath '/c3' fname ext]);
            i3 = spm_read_vols(HeaderInfo);

            mask = (i1+i2+i3)>0.5;
            mask = imdilate(mask,ones(settings.parameters.post.dilate,settings.parameters.post.dilate,settings.parameters.post.dilate));
            if settings.parameters.post.erode > 0
                mask = imerode(mask,ones(settings.parameters.post.erode,settings.parameters.post.erode,settings.parameters.post.erode));
            end

            HeaderInfo = spm_vol_nifti([fpath '/c1' fname ext]);

            HeaderInfo.fname = [fpath '/' fname '_brainmask' ext];
            HeaderInfo.private.dat.fname = HeaderInfo.fname;

            fprintf('Brain mask created.\n')
            spm_write_vol(HeaderInfo,mask);

            if ~settings.options.post.keepSegmentation
                delete([fpath '/c1' fname ext]);
                delete([fpath '/c2' fname ext]);
                delete([fpath '/c3' fname ext]);
            end
        end

        if settings.options.output.gzip
            gzip([fpath '/' fname ext]); delete([fpath '/' fname ext]);
            gzip([fpath '/' strrep([fname '_fwhm-' num2str(settings.parameters.post.fwhm) '_reg-' num2str(settings.parameters.post.reg) '_samp-' num2str(settings.parameters.post.samp)], '.', '') ext]); delete([fpath '/' strrep([fname '_fwhm-' num2str(settings.parameters.post.fwhm) '_reg-' num2str(settings.parameters.post.reg) '_samp-' num2str(settings.parameters.post.samp)], '.', '') ext]);
            if settings.options.post.mask && settings.options.post.keepSegmentation
                gzip([fpath '/' fname '_brainmask.nii']); delete([fpath '/' fname '_brainmask.nii']);
                gzip([fpath '/c1' fname ext]); delete([fpath '/c1' fname ext]);
                gzip([fpath '/c2' fname ext]); delete([fpath '/c2' fname ext]);
                gzip([fpath '/c3' fname ext]); delete([fpath '/c3' fname ext]);
            end
        end

        if settings.options.post.segmentation
            fprintf('| Biasfield correction and segmentation done.\n')
        else
            fprintf('| Biasfield correction done.\n')
        end
        fprintf('===============================================================\n');
    end
end
