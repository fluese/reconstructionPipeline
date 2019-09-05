function [bildcomp,phasemap]=tv_pc_silent(bildcomp,mu,runs,meanint)
% function [bildcomp,phasemap]=tv_pc_silent(bildcomp,mu,repetitions,intensity-rescaling)
%
% Function to extract the signal phase slice wise.
%
% Adjust mu to change phase map smoothness. (standard: mu = 0.1)
% Try values in the range of 0.05 to 0.2 for a start.

% Mu is set to work well with DICOM data. For simplicity, all data are
% temporarily rescaled to mimic DICOM scales [0 4095].

% lxm=round(size(bildcomp,1)/2);
% lym=round(size(bildcomp,2)/2);
% lzm=round(size(bildcomp,3)/2);
% if size(bildcomp,4)>1
%     lrm=round(size(bildcomp,4)/2);
% else
%     lrm=1;
% end

% bxo=squeeze(angle(bildcomp(lxm,:,:,lrm)));
% byo=squeeze(angle(bildcomp(:,lym,:,lrm)));
% bzo=squeeze(angle(bildcomp(:,:,lzm,lrm)));
% 
% axo=squeeze(abs(bildcomp(lxm,:,:,lrm)));
% ayo=squeeze(abs(bildcomp(:,lym,:,lrm)));
% azo=squeeze(abs(bildcomp(:,:,lzm,lrm)));

if nargin<2 || isempty(mu)
    mu=0.01;
end
if nargin<3 || isempty(runs)
    runs=1;
end
if nargin<4 || isempty(meanint)
    meanint=1000;
end

bildcomp=bildcomp.*meanint; % rescaling

% set up the deconvolution kernel
H       = fspecial('gaussian', [9 9], 0.1);

% Setup parameters (for example)
opts.rho_r   = 6;
opts.rho_o   = 100;
opts.beta    = [1 1 0];
opts.print   = false;
opts.alpha   = 0.7;
opts.method  = 'l2';

% Setup mu
% mu           = 0.1;

% Main routine
phasemap=zeros(size(bildcomp));
bildcomp_pha=zeros(size(bildcomp));
for lrun=1:runs
    for lr=1:size(bildcomp,4)
        bild=bildcomp(:,:,:,lr);
        bild_ang=zeros(size(bild));
        for lz=1:size(bild,3)
            g=bild(:,:,lz).*1;
            outr = deconvtv(real(g), H, mu, opts);
            outi = deconvtv(imag(g), H, mu, opts);
            bild_ang(:,:,lz)=angle(outr.f+1i.*outi.f);
        end
        bild_pc=abs(bild).*exp(1i.*(angle(bild)-bild_ang));
        bild_pc(bild==0)=0;
        bildcomp(:,:,:,lr)=single(bild_pc);
        bildcomp_pha(:,:,:,lr)=bild_ang;
    end
    phasemap=angle(exp(1i.*(phasemap+bildcomp_pha)));
end

bildcomp=bildcomp./meanint;

% bxp=squeeze(angle(bildcomp(lxm,:,:,lrm)));
% byp=squeeze(angle(bildcomp(:,lym,:,lrm)));
% bzp=squeeze(angle(bildcomp(:,:,lzm,lrm)));
% 
% bxe=squeeze(phasemap(lxm,:,:,lrm));
% bye=squeeze(phasemap(:,lym,:,lrm));
% bze=squeeze(phasemap(:,:,lzm,lrm));

% figure;
% subplot(3,1,1)
% imagesc(azo)
% axis image
% title('Magnitude of z-slice')
% subplot(3,1,2)
% imagesc([bzo bzp bze],[-pi pi])
% axis image
% colorbar
% title('Phase maps: Original --- Phase Corrected --- Extracted Background')
% subplot(3,1,3)
% rose(angle(bildcomp(:)),100);
% title('Histogram of phase over all volumes.')
% 
% disp('Please check the phase corrected object and the extracted phase map.')
% disp('The background in the extracted phase map should be smooth.')
% disp('The phase in the phase corrected object should be zero +- noise.')
% disp('If the background in the extracted phase is not smooth, try to decrease the intensity rescaling factor.')
% disp('If the phase of the phase corrected object is not properly corrected, try to increase the number of repetitions.')