clear all
close all
clc

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Demo file for deconvtv
% Image 'salt and pepper' noise removal
% 
% Stanley Chan
% University of California, San Diego
% 20 Jan, 2011
%
% Copyright 2011
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Prepare images
load lena
f_orig  = im2double(double(lena)./256);
[rows cols frames] = size(f_orig);
H       = fspecial('gaussian', [9 9], 0.1);
% g       = imfilter(f_orig, H, 'circular');
% g       = imnoise(g, 'salt & pepper', 0.05);

% Setup parameters (for example)
opts.rho_r   = 6;
opts.rho_o   = 100;
opts.beta    = [1 1 0];
opts.print   = true;
opts.alpha   = 0.7;
opts.method  = 'l2';

% Setup mu
mu           = 0.1;

% Main routine
tic
outr = deconvtv(real(g), H, mu, opts);
outi = deconvtv(imag(g), H, mu, opts);
out.f=outr.f+1i.*outi.f;
toc

bild=bildcomp(:,:,:,3);
bild_ang=zeros(size(bild));
for lz=1:size(bild,3)
    g=bild(:,:,lz).*10000;
    outr = deconvtv(real(g), H, mu, opts);
    outi = deconvtv(imag(g), H, mu, opts);
    bild_ang(:,:,lz)=angle(outr.f+1i.*outi.f);
end
bild_pc=abs(bild).*exp(1i.*(angle(bild)-bild_ang));
bild_pc(bild==0)=0;

% Display results
ccc=parula(512);
figure(1);
imagesx([angle(g) angle(out.f) angle(exp(1i.*(angle(g)-angle(out.f))))],[-pi pi]);
colormap(ccc);
title('input - output - difference');

% figure(2);
% subplot(1,2,2)
% imagesx(angle(out.f),[-pi pi]);
% colormap(ccc);
% title('output');