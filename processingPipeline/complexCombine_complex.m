close all; clear; clc;

% Change folder
cd D:\data\Reconstruction\yv98_4496\4Jahan

% Read in data
disp('Loading data.')
data1 = load_untouch_nii('yv98_4496-19_decorrelated-noiseScan_denoised-BM4D_distribution-Gauss_profile-manual_Tukey-1_unwrapped_maxChannel_adaptComb-10x10x10_complex.nii.gz');
data2 = load_untouch_nii('yv98_4496-21_denoised-BM4D_Tukey-1_complex_registerered_to_yv98_4496-19.nii.gz');
data3 = load_untouch_nii('yv98_4496-22_denoised-BM4D_Tukey-1_complex_registerered_to_yv98_4496-19.nii.gz');
data4 = load_untouch_nii('yv98_4496-23_denoised-BM4D_Tukey-1_complex_registerered_to_yv98_4496-19.nii.gz');

% Extract complex image information
disp('Extracting complex data.')
tmp1 = complex(data1.img(:,:,:,1),data1.img(:,:,:,2));
tmp2 = complex(data2.img(:,:,:,1),data2.img(:,:,:,2));
tmp3 = complex(data3.img(:,:,:,1),data3.img(:,:,:,2));
tmp4 = complex(data4.img(:,:,:,1),data4.img(:,:,:,2));

% Complex average
disp('Averaging data and creating complex, magntidude and phase data')
comp  = (tmp1+tmp2+tmp3+tmp4)/4;
mag   = abs((tmp1+tmp2+tmp3+tmp4)/4);
mag   = (mag/max(mag(:))).* 1000;
phase = angle((tmp1+tmp2+tmp3+tmp4)/4);

% Save as NIfTI file
disp('Saving data.')
data1=rmfield(data1,'untouch');
data1.img=comp;
save_nii(data1, 'yv98_4496-19+21+22+23_denoised-BM4D_Tukey-1_unwrapped_complex.nii.gz');

data1.hdr.dime.dim(5)=1;
data1.img=mag;
save_nii(data1, 'yv98_4496-19+21+22+23_denoised-BM4D_Tukey-1_unwrapped_magnitude.nii.gz');

data1.img=phase;
save_nii(data1, 'yv98_4496-19+21+22+23_denoised-BM4D_Tukey-1_unwrapped_phase.nii.gz');