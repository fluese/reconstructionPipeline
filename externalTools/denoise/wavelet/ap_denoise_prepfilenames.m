function para=ap_denoise_prepfilenames(para)


denames={'pkfile','denoise','wavecoeffs','parafile','phasemap'};



for la=1:numel(denames)
    para.(denames{la})=fullfile(para.storedir,[para.fileumbrella denames{la} '.mat']);
end
if(isfield(para,'refdat'))
    para.fileref=fullfile(para.storedir,para.refdat);
end


para.soso1=fullfile(para.storedir,[para.fileumbrella '_soso1' '.mat']);


para.sosd1=fullfile(para.storedir,[para.fileumbrella '_sosd1' '.mat']);


para.soso4=fullfile(para.storedir,[para.fileumbrella '_soso4' '.mat']);


para.sosd4=fullfile(para.storedir,[para.fileumbrella '_sosd4' '.mat']);