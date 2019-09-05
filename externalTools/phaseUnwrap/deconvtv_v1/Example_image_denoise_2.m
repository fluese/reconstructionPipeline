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
H       = fspecial('gaussian', [9 9], 1);
g       = imfilter(f_orig, H, 'circular');
% g       = imnoise(g, 'salt & pepper', 0.15);
g       = imnoise(g, 'localvar', 0.05.*f_orig);

% Setup parameters (for example)
opts.rho_r   = 5;
opts.rho_o   = 100;
opts.beta    = [1 1 0];
opts.print   = true;
opts.alpha   = 0.1;
opts.method  = 'l2';

% Setup mu
mu           = 10;

% Main routine
tic
out = deconvtv(g, H, mu, opts);
toc

% Display results
figure(1);
imshow(g);
title('input');

figure(2);
imshow(out.f);
title('output');