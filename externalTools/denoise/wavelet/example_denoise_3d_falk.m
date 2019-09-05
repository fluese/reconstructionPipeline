clear

%para.datfile='dat3d_5de';
% para.datfile='testdat';

para.storedir='D:\data\Reconstruction\yv98_4496\';
para.refdat='yv98_4496-19+21+22+23_average.nii.gz';

if(~exist('data','var'))
    disp('loading data ...')
    tstart=tic;
    bla=load_untouch_nii([para.storedir 'yv98_4496-19_decorrelated-noiseScan_Tukey-02_complex_perChannel.nii.gz']);
    % bla2=load_untouch_nii([para.storedir 'yv98_4496-19_decorrelated-noiseScan_Tukey-02_complex_perChannel.nii.gz']);
    % bla3=load_untouch_nii([para.storedir 'yv98_4496-19_decorrelated-noiseScan_Tukey-02_complex_perChannel.nii.gz']);
    % bla4=load_untouch_nii([para.storedir 'yv98_4496-19_decorrelated-noiseScan_Tukey-02_complex_perChannel.nii.gz']);
    tpk=toc(tstart);
    disp(['data loading duration :: ' datestr(tpk/(24*60*60), 'DD:HH:MM:SS.FFF')])

    data=zeros(size(bla.img,1),size(bla.img,2),size(bla.img,3),size(bla.img,4)/2,'single');
    for i=1:32
        data(:,:,:,i)=complex(bla.img(:,:,:,i),bla.img(:,:,:,32+i));
    end
%     global data;
%     data=tmp;

    clear bla i 
else
    disp('data are already in workspace!')
end

% para.fileumbrella='falk5';
% use a useful name
para.fileumbrella='nocheck';
para=ap_denoise_prepfilenames(para);

para.Ncoils=size(data,4); 
para.NRep=size(data,5);

% Definition of some parameters

% Phase correction parameters
para.pk.name={'gauss','gauss','gauss'};
para.pk.size=[size(data,1) size(data,2) size(data,3)];
para.pk.width=[8 8 8];
para.pk.shifts=[0 0 0];
para.pk.writepk=false;
para.pk.writephasemap=false;
para.pk.rescaled=true;
para.pk.scale=true;


para.pk.done=false; % This ensures pk is redone

tstart=tic;
[para, data]=ap_pk_3d(para, data);
tpk=toc(tstart);
disp(['Phase correction duration :: ' datestr(tpk/(24*60*60), 'DD:HH:MM:SS.FFF')])

para.shrinkcheck=false;

para.NwLevel=4;
para.wname='db5';

para.save.waves=true;
para.save.dwaves=true;
para.save.sos=true;
para.save.sosav=true;
para.save.nifti=true;

para=ap_denoise_3d(para, data);




%+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% WITH CHECKS
%+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

para.fileumbrella='withcheck';
para.shrinkcheck=true;
para.pk.done=true;

para=ap_denoise_prepfilenames(para);
para=ap_denoise_3d(para);



%AUSWERTUNG

datassets={'soso1','sosd1','soso4','sosd4'};

thisfile=para.fileref;

if(exist(thisfile,'file'))
    ref=load_untouch_nii(thisfile);
end
scf=max(ref.img(:));

mask=load_untouch_nii(fullfile(para.storedir,'brain.nii'));

para.fileumbrella='withcheck';
para=ap_denoise_prepfilenames(para);

clear withcheckdata xx
for la=1:numel(datassets)
    disp(para.(datassets{la}))
    xx=load(para.(datassets{la}));
    withcheckdata.(datassets{la})=scf*xx.(datassets{la})/max(xx.(datassets{la})(:));
    ref.imag=withcheckdata.(datassets{la});
    save_untouch_nii(ref,para.(datassets{la}))
end

para.fileumbrella='nocheck';
para=ap_denoise_prepfilenames(para);

clear nocheckdata xx
for la=1:numel(datassets)
    disp(para.(datassets{la}))
    xx=load(para.(datassets{la}));
    nocheckdata.(datassets{la})=scf*xx.(datassets{la})/max(xx.(datassets{la})(:));
    ref.imag=scf*nocheckdata.(datassets{la})/max(nocheckdata.(datassets{la})(:));
    save_untouch_nii(ref,para.(datassets{la}))
end

for la=1:numel(datassets)
nb.(datassets{la})=nocheckdata.(datassets{la});
wb.(datassets{la})=withcheckdata.(datassets{la});
end








