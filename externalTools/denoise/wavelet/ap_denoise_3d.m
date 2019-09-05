function para=ap_denoise_3d(para, data)

% global data
% global dataDEN
% 
% if(isempty(data))
%     error('THE GLOBAL VARIABLE data IS NOT DEFINED!!!')
% end




[wpc, Nlevel]= ap_wavedec3d(para, data);
para         = ap_wavenoise(para, wpc,Nlevel);
wpc          = ap_waveden(para, wpc, Nlevel, size(data,4));
dataDEN      = ap_waverec3d(para, size(data), wpc);
ap_save(para, data, dataDEN);


    function [wpc,Nlevel]=ap_wavedec3d(para, data)
        wpc=cell(para.Ncoils,para.NRep);
        tstart=tic;
        for icoil=1:para.Ncoils
            for irep=1:para.NRep
                wpc{icoil,irep} = wavedec3(data(:,:,:,icoil,irep),para.NwLevel,para.wname,'zpd');
            end
        end
        twd=toc(tstart);
        disp(['3d WaveDec duration :: ' datestr(twd/(24*60*60), 'DD:HH:MM:SS.FFF')])
        Nlevel=numel(wpc{1,1}.dec);
    end

    function para=ap_wavenoise(para, wpc,Nlevel)
        wnoisestdvecre=zeros(para.Ncoils,para.NRep,Nlevel);
        wnoisestdvecim=zeros(para.Ncoils,para.NRep,Nlevel);
        
        tstart=tic;
        for icoil=1:para.Ncoils
            for irep=1:para.NRep
                for ilevel=1:Nlevel
                    wnoisestdvecre(icoil,irep,ilevel)=mad(real(wpc{icoil,irep}.dec{ilevel}(:)),1)*1.4826;
                    wnoisestdvecim(icoil,irep,ilevel)=mad(imag(wpc{icoil,irep}.dec{ilevel}(:)),1)*1.4826;
                end
            end
        end
        para.wnoisestd=min([wnoisestdvecim(:);wnoisestdvecre(:)]);
        twd=toc(tstart);
        disp(['Noise estimation duration :: ' datestr(twd/(24*60*60), 'DD:HH:MM:SS.FFF')])
    end

    function wpc=ap_waveden(para, wpc, Nlevel, sz)
        wnoisestd=para.wnoisestd;
        shrinkcheck=para.shrinkcheck;
        %wpcd=wpc;
        % Maybe later in case wpc and wpcd shall be analzyed
        % separately
        for irep=1:para.NRep
            disp(['NREP=' num2str(irep)])
            tstart=tic;
            for ilevel=2:Nlevel
                aa=zeros([size(wpc{1,1}.dec{ilevel}) sz], 'single');
                
                for icoil=1:size(wpc,1)
                    aa(:,:,:,icoil)=wpc{icoil,irep}.dec{ilevel};
                end
                
                bb=zeros(size(aa), 'single');
                
                for ix=1:size(aa,1)
                    
                    for iy=1:size(aa,2)
                        for iz=1:size(aa,3)
                            vec=aa(ix,iy,iz,:);
                            
                            rdenvec=stein(real(vec(:)),wnoisestd);
                            idenvec=1i*stein(imag(vec(:)),wnoisestd);
                            if(shrinkcheck)
                                checkdenreal=std(rdenvec)/std(real(vec));
                                checkdenimag=std(idenvec)/std(imag(vec));
                                if(checkdenreal>0.7)
                                    rdenvec=real(vec(:));
                                end
                                if(checkdenimag>0.7)
                                    idenvec=imag(vec(:));
                                end
                            end
                            denvec=rdenvec+1i*idenvec;
                            
                            bb(ix,iy,iz,:)=denvec;
                            
                        end
                    end
                end
                
                for icoil=1:size(wpc,1)
                    wpc{icoil,irep}.dec{ilevel}=bb(:,:,:,icoil);
                end
                
            end
            twd=toc(tstart);
            disp(['Denosing duration :: ' datestr(twd/(24*60*60), 'DD:HH:MM:SS.FFF')])
        end
        
        % TOO SLOOOOOOOOOOW!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        
%         for irep=1:para.NRep
%             disp(['NREP=' num2str(irep)])
%             
%             tic
%             for ilevel=2:Nlevel
%                 
%                 SIx=size(wpc{1,1}.dec{ilevel},1);
%                 SIy=size(wpc{1,1}.dec{ilevel},2);
%                 SIz=size(wpc{1,1}.dec{ilevel},3);
%                 
%                 for ix=1:SIx
%                     for iy=1:SIy
%                         for iz=1:SIz
%                             
%                             vec=zeros(size(wpc,irep),1);
%                             for ichannel=1:size(wpc,1)
%                                 vec(ichannel)=wpc{ichannel,irep}.dec{ilevel}(ix,iy,iz);
%                             end
%                             shrinkto=zeros(size(vec));
%                             rdenvec=stein(real(vec),para.wnoisestd,shrinkto);
%                             idenvec=1i*stein(imag(vec),para.wnoisestd,shrinkto);
%                             if(para.shrinkcheck)
%                                 checkdenreal=std(rdenvec)/std(real(vec));
%                                 checkdenimag=std(idenvec)/std(imag(vec));
%                                 if(checkdenreal>0.7)
%                                     rdenvec=real(vec);
%                                 end
%                                 if(checkdenimag>0.7)
%                                     idenvec=imag(vec);
%                                 end
%                             end
%                             denvec=rdenvec+1i*idenvec;
%                             
%                             
%                             for ichannel=1:size(wpc,1)
%                                 wpc{ichannel,irep}.dec{ilevel}(ix,iy,iz)=denvec(ichannel);
%                             end
%                             
%                         end
%                     end
%                 end
%                 
%             end
%             toc
%         end
        
    end

    function dataDEN=ap_waverec3d(para, sz, wpc)
        dataDEN=complex(zeros(sz,'single'));
        
        tstart=tic;
        for icoil=1:para.Ncoils
            for irep=1:para.NRep
                dataDEN(:,:,:,icoil,irep)=waverec3(wpc{icoil,irep});
            end
            
        end
        twd=toc(tstart);
        disp(['3D WaveRec duration :: ' datestr(twd/(24*60*60), 'DD:HH:MM:SS.FFF')])
    end

    function ap_save(para, data, dataDEN)
        if(para.save.nifti)
            bildosos = squeeze(sum(abs(data.^2),4)).^(1/2);
            bilddsos = squeeze(sum(abs(dataDEN.^2),4)).^(1/2);
            bildo4avcomp=sum(data,5)/4;
            bildd4avcomp=sum(dataDEN,5)/4;
            
            bildo4avcompsos=squeeze(sum(abs(bildo4avcomp.^2),4)).^(1/2);
            bildd4avcompsos=squeeze(sum(abs(bildd4avcomp.^2),4)).^(1/2);
            
            soso1=single(1e20*bildosos(:,:,:,1)/max(bildosos(:)));
            sosd1=single(1e20*bilddsos(:,:,:,1)/max(bilddsos(:)));
            soso4=single(1e20*bildo4avcompsos/max(bildo4avcompsos(:)));
            sosd4=single(1e20*bildd4avcompsos/max(bildd4avcompsos(:)));
            
            save(para.soso1,'soso1')
            save(para.sosd1,'sosd1')
            save(para.soso4,'soso4')
            save(para.sosd4,'sosd4')
            
            %             thisfile=para.fileref;
            %             if(exist(thisfile,'file'))
            %                 ref=load_untouch_nii(thisfile);
            %             end
            %             save_untouch_nii(ref,para.reffile)
            %
            %             % single channel
            %             soso1=ref;
            %             soso1.img=single(1000*bildosos(:,:,:,1)/max(bildosos(:)));
            %             save_untouch_nii(soso1,para.soso1)
            %
            %             sosd1=ref;
            %             sosd1.img=single(1000*bilddsos(:,:,:,1)/max(bilddsos(:)));
            %             save_untouch_nii(sosd1,para.sosd1)
            %
            %             % 4 NAV
            %             soso4=ref;
            %             soso4.img=single(1000*bildo4avcompsos/max(bildo4avcompsos(:)));
            %             save_untouch_nii(soso4,para.soso4)
            %
            %
            %             sosd4=ref;
            %             sosd4.img=single(1000*bildd4avcompsos/max(bildd4avcompsos(:)));
            %             save_untouch_nii(sosd4,para.sosd4)
        end
    end

end