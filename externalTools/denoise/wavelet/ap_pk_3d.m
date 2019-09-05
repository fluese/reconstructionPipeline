function [para,data,phasemap]=ap_pk_3d(para,data)
% Phase corecction of 3D data
% Phasemaps are derived from
% windowed data and wavelet thresholding.
%
if(para.pk.done)
    warning('Nothing done. data were already phasecorrected!')
    phasemap=nan;
    return
end
% global data
% if(isempty(data))
%     error('THE GLOBAL VARIABLE data IS NOT DEFINED!!!')
% end

Ncoils=size(data,4); NRep=size(data,5);

wt=single(mrt_wdw(para.pk));

bk=complex(zeros(size(data),'single'));

for icoil=1:Ncoils
    for irep=1:NRep
        bk(:,:,:,icoil,irep)=ifftshift(ifftshift(ifftshift(ifft(ifft(ifft(fftshift(fftshift(fftshift(data(:,:,:,icoil,irep),1),2),3),[],1),[],2),[],3),1),2),3).*wt;
        bk(:,:,:,icoil,irep)=fftshift(fftshift(fftshift(fft(fft(fft(ifftshift(ifftshift(ifftshift(bk(:,:,:,icoil,irep),1),2),3),[],1),[],2),[],3),1),2),3);
    end
end

clear wt
ws=cell(Ncoils,NRep);

for icoil=1:Ncoils
    for irep=1:NRep
        ws{icoil,irep} = wavedec3(bk(:,:,:,icoil,irep),4,'db5','zpd');
        sz=numel(ws{1,1}.dec(:));
        for i=1:sz
            ws{icoil,1}.dec{i}=single(ws{icoil,1}.dec{i});
        end
    end
end

wsclean=ws;
Nlevel=numel(ws{1,1}.dec);

%     sbk=complex(zeros(size(bk)),'single');


for icoil=1:Ncoils
    for irep=1:NRep
        %         im_W=ws{icoil,irep}.dec{29};
        %         m = sort(abs(im_W(:)),'descend');
        %         % Only 10% are kept!
        %         % nice way to compute threshold
        %         ndx = floor(length(m)*0.01/100);
        %         thresh = m(ndx)
        %         im_W_th  = im_W .* (abs(im_W) > thresh);
        for ilevel=8:Nlevel
            wsclean{icoil,irep}.dec{ilevel}=ws{icoil,irep}.dec{ilevel}*0;
        end
        bk(:,:,:,icoil,irep)=waverec3(wsclean{icoil,irep});
    end

end

clear ws wsclean

phasemap=angle(bk);
clear bk

% if(~exist('data','var'))
%     disp('loading data ...')
%     tstart=tic;
%     bla=load_untouch_nii([para.storedir 'yv98_4496-19_decorrelated-noiseScan_Tukey-02_complex_perChannel.nii.gz']);
%     % bla2=load_untouch_nii([para.storedir 'yv98_4496-19_decorrelated-noiseScan_Tukey-02_complex_perChannel.nii.gz']);
%     % bla3=load_untouch_nii([para.storedir 'yv98_4496-19_decorrelated-noiseScan_Tukey-02_complex_perChannel.nii.gz']);
%     % bla4=load_untouch_nii([para.storedir 'yv98_4496-19_decorrelated-noiseScan_Tukey-02_complex_perChannel.nii.gz']);
%     tpk=toc(tstart);
%     disp(['data loading duration :: ' datestr(tpk/(24*60*60), 'DD:HH:MM:SS.FFF')])
% 
%     data=zeros(size(bla.img,1),size(bla.img,2),size(bla.img,3),size(bla.img,4)/2,'single');
%     for i=1:32
%         data(:,:,:,i)=complex(bla.img(:,:,:,i),bla.img(:,:,:,32+i));
%     end
% %     global data;
% %     data=tmp;
% 
%     clear bla i 
% else
%     disp('data are already in workspace!')
% end

% global data
if(para.pk.scale)
    scalefac=max(abs(data(:)));
    data=(4096/scalefac)*data.*exp(-1i*phasemap);
    para.pk.scalefac=scalefac;
    para.pk.rescaled=true;
else
    data=data.*exp(-1i*phasemap);
    para.pk.rescaled=false;
end
para.pk.done=true;

if(para.pk.writepk)
    save(para.pkfile,'data','para','-v7.3');
end
if(para.pk.writephasemap)
    save(para.phasemap,'phasemap','para','-v7.3');
end


