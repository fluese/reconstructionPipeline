% Function to run a simple quality assessment on reconstructed data.
%
% ************************************************************************
% Version 1.0                                                  15.02.2019
% Falk Luesebrink              falk dot luesebrink at med dot ovgu dot de
% ************************************************************************
clear
clc
close all

[fileRef, pathRef] = uigetfile('*.nii;*.nii.gz', 'Select reference file.', 'MultiSelect', 'off');
[files, path] = uigetfile('*.nii;*.nii.gz', 'Select files which shall be quality checked.', 'MultiSelect', 'on');

ref=load_untouch_nii([pathRef fileRef]);
ref.img = uint8(ref.img);
% ref.img=uint8((ref.img./max(ref.img(:))).*255);

% useNoise = true;
% if useNoise
%     userdir = 'D:\MR\';
%     if ~exist('userdir', 'var')
%         if ispc
%           userdir = getenv('USERPROFILE');
%         else
%           userdir = getenv('HOME');
%         end
%     end
%     pathNoise = uigetdir(userdir,'Select folder were noise data is stored in.');
%     DIRNAMES  = dir(pathNoise);
% 
%     sz = size(ref.img);
%     noiseData = complex(zeros(sz(1), sz(2), length(DIRNAMES)-2, 'single'));
%     for Channel = 1:length(DIRNAMES)-2
%         tmp = load([pathNoise '/Channel_' num2str(Channel) '.mat']);
%         noiseData(:,:,length(DIRNAMES)-2) = tmp.data;
%     end
%     noiseData = squeeze(sum(abs(noiseData.^2),3)).^(1/2);
%     % noise.img=uint8((noise.img./max(noise.img(:))).*255);
% end

if ~iscell(files)
    tmp=files;
    files=cell(1,1);
    files{1,1}=tmp;
    clear tmp
end

PSNR       = zeros(1,length(files));
% SNR        = zeros(1,length(files));
SSIM       = zeros(1,length(files));

for m = 1:length(files)
    fprintf('Processing file: %s\n', files{1,m});
    data=load_untouch_nii([path files{1,m}]);
%     data.img=uint8((data.img./max(data.img(:))).*255);
    data.img = uint8(data.img);
    % Calculate the metrics
%     fprintf('Calculating SNR...\n');
%     signal = zeros(25,30, 'single');
%     noise  = zeros(25,30, 'single');
%     for x = 201:225
%         for y = 201:230
%             signal(x-200,y-200) = data.img(x,y,160);
%             noise(x-200,y-200) = noiseData(x,y);
%         end
%     end
%     SNR(m) = mean(signal(:))/std(noise(:));
    fprintf('Calculating PSNR...\n');
    PSNR(m) = psnr(data.img(:,:,90),ref.img(:,:,90));
%     PSNR(m) = psnr(data.img,ref.img);
    fprintf('Calculating SSIM...\n');
    SSIM(m) = ssim(data.img(:,:,90),ref.img(:,:,90));  
%     SSIM(m) = ssim(data.img,ref.img); 
%     fprintf('Calculating ESP...\n');
%     hold on
%     esp(ref.img(:,:,140), data.img(:,:,140), 0.7);
    fprintf('Done\n\n');
end
% hold off

%     %% Save metrics in files
%       fname = sprintf('PSNR(d= %d).mat', d);
%     save(fname, 'PSNR');
%       fname = sprintf('PSNRHVSM(d= %d).mat',d);
%     save(fname, 'PSNRHVSM');
%       fname = sprintf('SSIM(d= %d).mat', d);
%     save(fname, 'SSIM');
%       fname = sprintf('SNR(d= %d).mat', d);
%     save(fname, 'SNR');
    
    %% plot

% %Plot PSNR
% figure
% plot (M_vec,PSNR,'-s')
% legend(sprintf('Patchgroe�e = %g', d)) 
% title('Influence of the search radius on PSNR')
% xlabel('search radius') % x-axis label
% ylabel('PSNR') % y-axis label
% grid on
% %Plot MSSIM
% figure
% 
% plot (M_vec,SSIM,'-s')
% legend(sprintf('Patchgroe�e = %g', d)) 
% title('Influence of the search radius on MSSIM and SSIM')
% xlabel('search radius') % x-axis label
% ylabel('MSSIM') % y-axis label
% grid on
% %Plot PSNRHVSM
% figure
% plot (M_vec,PSNRHVSM,'-s')
% legend(sprintf('Patchgroe�e = %g', d)) 
% title('Influence of the search radius on PSNRHVSM')
% xlabel('search radius') % x-axis label
% ylabel('PSNRHVSM') % y-axis label
% grid on
% 
% %Plot SNR
% figure
% plot (M_vec,SNR,'-s')
% legend(sprintf('Patchgroe�e = %g', d)) 
% title('Influence of the search radius on SNR')
% xlabel('search radius') % x-axis label
% ylabel('SNR') % y-axis label
% grid on