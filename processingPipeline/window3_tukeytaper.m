% This function creates a 2 dimentional window for a sample image, it takes
% the dimension of the window and applies the 1D window function
% This is does NOT using a rotational symmetric method to generate a 2 window
%
% Disi A ---- May,16, 2013
%
% Modified by Falk Luesebrink to allow creation of a 3D window and tapering of
% Tukey window. (19.04.2018) - falk dot luesebrink at med dot ovgu dot de
%
%     [N,M]=size(imgage);
% ---------------------------------------------------------------------
%     w_type is defined by the following 
%     @bartlett       - Bartlett window.
%     @barthannwin    - Modified Bartlett-Hanning window. 
%     @blackman       - Blackman window.
%     @blackmanharris - Minimum 4-term Blackman-Harris window.
%     @bohmanwin      - Bohman window.
%     @chebwin        - Chebyshev window.
%     @flattopwin     - Flat Top window.
%     @gausswin       - Gaussian window.
%     @hamming        - Hamming window.
%     @hann           - Hann window.
%     @kaiser         - Kaiser window.
%     @nuttallwin     - Nuttall defined minimum 4-term Blackman-Harris window.
%     @parzenwin      - Parzen (de la Valle-Poussin) window.
%     @rectwin        - Rectangular window.
%     @taylorwin      - Taylor window.
%     @tukeywin       - Tukey window.
%     @triang         - Triangular window.
%
%   Example: 
%   To compute windowed 2D fFT
%   [r,c]=size(img);
%   w=window2(r,c,@hamming);
% 	fft2(img.*w);

function w=window3_tukeytaper(sz,w_func,taper)

N=sz(1);
M=sz(2);
L=sz(3);

wc=window(w_func,N,taper);
wr=window(w_func,M,taper);
ws=window(w_func,L,taper);
[maskr,maskc,masks]=meshgrid(wr,wc,ws);

clear wc wr ws

w=single(maskr.*maskc.*masks);
end

