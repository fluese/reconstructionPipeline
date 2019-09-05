function data=phaseCorrection(name, settings, data, Channel)
% Phase corecction of 3D data
% Phasemaps are derived from windowed data and wavelet thresholding.
%
% By Andre Pampel (MPI Leipzig)

sz=[size(data,1) size(data,2) size(data,3)];
wt=mrt_wdw(settings.parameters.pk, sz);
    
bk=fftshift(fftn(ifftshift(data.*wt)));

% figure; imagesc(angle(squeeze(bk(:,:,20))));

ws = wavedec3(real(bk),4,'db5','zpd');
Nlevel=numel(ws.dec);

%         im_W=ws{icoil,irep}.dec{29};
%         m = sort(abs(im_W(:)),'descend');
%         % Only 10% are kept!
%         % nice way to compute threshold
%         ndx = floor(length(m)*0.01/100);
%         thresh = m(ndx)
%         im_W_th  = im_W .* (abs(im_W) > thresh);
for ilevel=8:Nlevel
    ws.dec{ilevel}=ws.dec{ilevel}*0;
end
bk=waverec3(ws);

clear ws

phasemap=angle(bk);

figure; imagesc(angle(squeeze(phasemap(:,:,20))));

% % Phase correction parameters
% 
% 
% phasemap=data.*wt;
% clear wt
% phasemap=fftshift(fftshift(fftshift(fft(fft(fft(ifftshift(ifftshift(ifftshift(phasemap,1),2),3),[],1),[],2),[],3),1),2),3);
% 
% ws = wavedec3(phasemap,4,'db5','zpd');
% 
% wsclean=ws;
% Nlevel=numel(ws.dec);
% 
% %         im_W=ws{icoil,irep}.dec{29};
% %         m = sort(abs(im_W(:)),'descend');
% %         % Only 10% are kept!
% %         % nice way to compute threshold
% %         ndx = floor(length(m)*0.01/100);
% %         thresh = m(ndx)
% %         im_W_th  = im_W .* (abs(im_W) > thresh);
% for ilevel=8:Nlevel
%     wsclean.dec{ilevel}=ws.dec{ilevel}*0;
% end
% phasemap=single(waverec3(wsclean));
% clear ws wsclean
% 
% phasemap=angle(phasemap);

if(settings.parameters.pk.scale)
    data=fftshift(fftshift(fftshift(fft(fft(fft(ifftshift(ifftshift(ifftshift(data,1),2),3),[],1),[],2),[],3),1),2),3);
    scalefac=max(abs(data(:)));
    data=(4096/scalefac)*data.*exp(-1i*phasemap);
    data=ifftshift(ifftshift(ifftshift(ifft(ifft(ifft(fftshift(fftshift(fftshift(data,1),2),3),[],1),[],2),[],3),1),2),3);
else
    data=fftshift(fftshift(fftshift(fft(fft(fft(ifftshift(ifftshift(ifftshift(data,1),2),3),[],1),[],2),[],3),1),2),3);
    data=data.*exp(-1i*phasemap);
    data=ifftshift(ifftshift(ifftshift(ifft(ifft(ifft(fftshift(fftshift(fftshift(data,1),2),3),[],1),[],2),[],3),1),2),3);
end

if(settings.parameters.pk.writephasemap)
    save([name.path 'phasemap' num2str(Channel) '.mat'],'phasemap','-v7.3');
end