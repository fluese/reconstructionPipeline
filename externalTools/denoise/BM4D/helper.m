classdef helper
    methods (Static)
        
        function [ map ] = getNoiseMap( y, noise_factor)
            % noise modulation field
            map = ones(3,3,3);
            map(2,2,2) = noise_factor;
            s = size(y);
            [x1,y1,z1] = meshgrid(1:3,1:3,1:3);
            [x2,y2,z2] = meshgrid(1:2/(s(2)-1):3,1:2/(s(1)-1):3,1:2/(s(3)-1):3);
            map = interp3(x1,y1,z1,map,x2,y2,z2,'cubic');
        end
        
        function [ C ] = constantsSparseTraj3D
            % constants
            C.RADIAL     = 1;
            C.SPIRAL     = 2;
            C.LOG_SPIRAL = 3;
            C.LIM_ANGLE  = 4;
            C.SPHERICAL  = 5;
            C.COMPLETE   = 6;
            C.RAND_LINES = 7;
            C.RAND_PTS   = 8;
            C.LOW_PASS   = 9;
            C.HELICAL    = 10;
            C.TRAJ = {'RADIAL', 'SPIRAL', 'LOG_SPIRAL', 'LIM_ANGLE', 'SPHERICAL', 'COMPLETE'};

            C.BRAINWEB   = 1;
            C.SHEPPLOGAN = 2;
            C.DATA = {'BRAINWEB','SHEPP-LOGAN'};

            C.REAL    = 1;
            C.COMPLEX = 2;
            C.OBS     = {'REAL','COMPLEX'};

            C.NONE  = 0;
            C.TEXT  = 1;
            C.IMAGE = 2;

            C.REC_2D = 0;
            C.REC_3D = 1;
        end
        
        function visualizeEstMap( y, sigma_est, eta )
            figure,
            subplot(1,2,1)
            hold all
            hist((sigma_est(:)-eta(:))./eta(:),40)
            h = findobj(gca,'Type','patch');
            set(h,'FaceColor','w','EdgeColor','k')
            line([0 0],get(gca,'ylim'), 'Color','k','LineWidth',2);
            box
            axis square
            ylabel('Voxels','Interpreter','Latex')
            title('Relative noise estimation error')
            hold off
            f = round(size(y,3)/2);
            subplot(1,2,2)
            img = y(:,:,f);
            map = abs((sigma_est(:,:,f)-eta(:,:,f))./eta(:,:,f));
            mask = cat(3, ones(size(img)), zeros(size(img)), zeros(size(img)));
            imshow(img,[],'InitialMagnification','fit')
            hold all
            h = imshow(mask,[]);
            set(h, 'AlphaData', map)
            title(['Estimation error map / cross-section ',int2str(f)])
            hold off
        end
        
        function visualizeXsect( varargin )
            if ~exist('crop','var')
                crop = 1;
            end
            
            if ndims(varargin{1})~=3
                error('Input must be 3-D.')
            end
            
            labels = cell(1,nargin);
            for l=1:length(labels)
                labels{l} = strrep(inputname(l),'_','\_');
            end
            
            % show 2-D cross-sections
            figure,
            gray = mean(get(gcf,'color'));
            f = round(size(varargin{1})/2);
            orientation = {'Horizontal','Coronal','Sagittal'};
            for d=1:ndims(varargin{1})
                subplot(3,1,d)
                dim = circshift(1:ndims(varargin{1}),[0,d-1]);
                
                data = permute(varargin{1}, dim);
                s1 = size(data,1+mod(ndims(varargin{1})-d+1,2));
                s2 = size(data,2-mod(ndims(varargin{1})-d+1,2));
                sb = round(max(size(varargin{1}))/6);
                band = gray*ones(s1,sb);
                
                img = [];
                for n=1:nargin
                    data = permute(varargin{n}, dim);
                    img = [img rot90(data(:,:,f(dim(3))),ndims(varargin{n})-d+1)];
                    if n<nargin
                        img = [img band];
                    end
                end
                imshow(img, [0 1], 'InitialMagnification','fit')
                title([orientation{d},' view / slice ',int2str(f(dim(3)))], ...
                    'Interpreter','latex')
                for l=1:length(labels)
                    text((l-1)*(s2+sb)+s2/2,s1, labels{l}, ...
                        'VerticalAlignment','top','HorizontalAlignment','center')
                end
            end
            
            % show 3-D cross-sections
            h = figure;
            %set(h, 'color', 'white');
            set(h, 'renderer', 'zbuffer');
            for i=1:length(varargin)
                y = varargin{i};
                if isempty(y)
                    continue;
                end
                x = round(size(y)./2);
                
                p1 = squeeze(y(x(1),:,:));
                p2 = squeeze(y(:,x(2),:))';
                p3 = squeeze(y(:,:,x(3)));
                
                if crop
                    offset = 1;
                    p1(:,1:x(3)-offset) = [];
                    p1(1:x(2)-offset,:) = [];
                    
                    p3(:,1:x(2)-offset) = [];
                    p3(1:x(1)-offset,:) = [];
                    
                    p2(:,1:x(1)-offset) = [];
                    p2(1:x(3)-offset,:) = [];
                    
                    x(:) = offset;
                end
                
                x1 = zeros(size(p1)) + x(1);
                x2 = zeros(size(p2)) + x(2);
                x3 = zeros(size(p3)) + x(3);
                [X1_1,X1_2] = meshgrid(1:size(p1,2),1:size(p1,1));
                [X2_1,X2_2] = meshgrid(1:size(p2,2),1:size(p2,1));
                [X3_1,X3_2] = meshgrid(1:size(p3,2),1:size(p3,1));
                
                subplot(1,length(varargin),i);
                hold all
                warp(X1_2,x1,X1_1,p1)
                warp(x2,X2_1,X2_2,p2)
                warp(X3_1,X3_2,x3,p3)
                view([35 30])
                axis off vis3d;
                camproj('persp')
                title(labels{i},'fontsize',12);
            end
        end
        
        function [ S ] = sampling( trajectory, n, rot_deg, line_nbr, line_std, coverage, tol )
            
            % load constants
            C = self_constantsSparseTraj3D;
            
            phantom_size = [n n n];
            S = ones(phantom_size);
            
            if trajectory==C.RADIAL || trajectory==C.LIM_ANGLE
                history = zeros(2);
                stop = 0;
                % getting line-mask (radial lines)
                S = zeros(phantom_size);
                scale = n*0.2;
                % scaling parameter to the desired coverage
                while (sum(sum(S(:,:,1)))==0 || abs(sum(sum(S(:,:,1)))/numel(S(:,:,1))-coverage)>tol) && ~stop
                    if sum(sum(S(:,:,1)))/numel(S(:,:,1))>coverage
                        scale = scale-tol/2;
                    else
                        scale = scale+tol/2;
                    end
                    if scale==history(1,1)
                        scale = history(1,history(2,:)==min(history(2,:)));
                        stop = 1;
                    end
                    S(:,:,1) = ifftshift( radial( trajectory, n, scale, 1, rot_deg, line_nbr, line_std ) );
                    history = circshift(history,[0 1]);
                    history(1,end) = scale;
                    history(2,end) = abs(sum(sum(S(:,:,1)))/numel(S(:,:,1))-coverage);
                end
                for i2=2:n
                    S(:,:,i2) = ifftshift( radial( trajectory, n, scale, i2, rot_deg, line_nbr, line_std ) );
                end
            end
            
            if trajectory==C.SPIRAL
                history = zeros(2);
                stop = 0;
                % getting line-mask (spiral lines)
                S = zeros(phantom_size);
                scale = n*6e-3;
                t = linspace(0,10*n,5e5);
                % scaling parameter to the desired coverage
                while (sum(sum(S(:,:,1)))==0 || abs(sum(sum(S(:,:,1)))/numel(S(:,:,1))-coverage)>tol) && ~stop
                    if sum(sum(S(:,:,1)))/numel(S(:,:,1))>coverage
                        scale = scale+1e-2;
                    else
                        scale = scale-1e-2;
                    end
                    if scale==history(1,1)
                        scale = history(1,history(2,:)==min(history(2,:)));
                        stop = 1;
                    end
                    S(:,:,1) = ifftshift( spiral( t, n, scale, 1, rot_deg, line_nbr, line_std ));
                    history = circshift(history,[0 1]);
                    history(1,end) = scale;
                    history(2,end) = abs(sum(sum(S(:,:,1)))/numel(S(:,:,1))-coverage);
                end
                for i2=2:n
                    S(:,:,i2) = ifftshift(spiral( t, n, scale, i2, rot_deg, line_nbr, line_std ));
                end
            end
            
            if trajectory==C.LOG_SPIRAL
                history = zeros(2);
                stop = 0;
                % getting line-mask (log-spiral lines)
                S = zeros(phantom_size);
                scale = 3/n;
                t = linspace(0,10*n,5e5);
                % scaling parameter to the desired coverage
                while (sum(sum(S(:,:,1)))==0 || abs(sum(sum(S(:,:,1)))/numel(S(:,:,1))-coverage)>tol) && ~stop
                    if sum(sum(S(:,:,1)))/numel(S(:,:,1))>coverage
                        scale = scale+1e-3;
                    else
                        scale = scale-1e-3;
                    end
                    if scale==history(1,1)
                        scale = history(1,history(2,:)==min(history(2,:)));
                        stop = 1;
                    end
                    S(:,:,1) = ifftshift( log_spiral( t, n, scale, 1, rot_deg, line_nbr, line_std ));
                    history = circshift(history,[0 1]);
                    history(1,end) = scale;
                    history(2,end) = abs(sum(sum(S(:,:,1)))/numel(S(:,:,1))-coverage);
                end
                for i2=2:n
                    S(:,:,i2) = ifftshift( log_spiral( t, n, scale, i2, rot_deg, line_nbr, line_std ) );
                end
            end
            
            if trajectory==C.HELICAL
                scale = 1;
                t = 0:n*1000;
                r = n/2;% 0.01*t;
                x = r.*cos(scale*2*pi*t/10000);
                y = r.*sin(scale*2*pi*t/10000);
                z = t/1000;
                
                temp = zeros([n n n]);
                S = zeros([n n n]);
                for i=1:length(t)
                    xr = max(1,min(n,round(x(i)+n/2)));
                    yr = max(1,min(n,round(y(i)+n/2)));
                    zr = max(1,min(n,round(z(i))));
                    
                    if temp(yr,xr,zr)==0
                        temp(yr,xr,zr) = 1;
                        xr = max(1,min(n,round(xr+line_std*randn)));
                        yr = max(1,min(n,round(yr+line_std*randn)));
                        zr = max(1,min(n,round(zr+line_std*randn)));
                        S(yr,xr,zr) = 1;
                    end
                end
                for i=1:size(S,3)
                    S(:,:,i) = ifftshift(S(:,:,i));
                end
                for i=1:size(S,1)
                    for j=1:size(S,2)
                        S(i,j,:) = squeeze(ifftshift(S(i,j,:)));
                    end
                end
            end
            
            if trajectory==C.SPHERICAL
                history = zeros(2);
                stop = 0;
                S = zeros(phantom_size);
                scale = 8;
                p0 = [0 0 0];
                p1 = size(S);
                t = 0:1/n:1;
                while (sum(S(:))==0 || abs(mean(S(:))-coverage)>tol) && ~stop
                    if mean(S(:))>coverage
                        scale = scale+1;
                    else
                        scale = scale-1;
                    end
                    if scale==history(1,1)
                        scale = history(1,history(2,:)==min(history(2,:)));
                        stop = 1;
                    end
                    S = zeros(phantom_size);
                    for sx=0:scale:p1(2)
                        for sy=0:scale:p1(1)
                            for sz=0:scale:p1(3)
                                a = p0;
                                b = p1;
                                a(1) = max(1,min(p1(2),a(1)+sx));
                                a(2) = max(1,min(p1(1),a(2)+sy));
                                a(3) = max(1,min(p1(3),a(3)+sz));
                                
                                b(1) = max(1,min(p1(2),b(1)-sx));
                                b(2) = max(1,min(p1(1),b(2)-sy));
                                b(3) = max(1,min(p1(3),b(3)-sz));
                                d = b-a;
                                x = max(1,min(p1(2),round(a(1)+d(1).*t+line_std*randn)));
                                y = max(1,min(p1(1),round(a(2)+d(2).*t+line_std*randn)));
                                z = max(1,min(p1(3),round(a(3)+d(3).*t+line_std*randn)));
                                for i=1:length(t)
                                    S(y(i),x(i),z(i)) = 1;
                                end
                            end
                        end
                    end
                    history = circshift(history,[0 1]);
                    history(1,end) = scale;
                    history(2,end) = abs(mean(S(:))-coverage);
                end
                for i=1:size(S,3)
                    S(:,:,i) = ifftshift(S(:,:,i));
                end
                for i=1:size(S,1)
                    for j=1:size(S,2)
                        S(i,j,:) = squeeze(ifftshift(S(i,j,:)));
                    end
                end
            end
            
            if trajectory==C.RAND_LINES
                S = zeros(phantom_size);
                p = coverage;
                for i2=1:n
                    l = rand(1,n);
                    S(l<p,:,i2) = 1;
                    S(:,:,i2) = ifftshift(S(:,:,i2));
                end
            end
            
            if trajectory==C.RAND_PTS
                S = zeros(phantom_size);
                p = coverage;
                for i2=1:n
                    l = rand(n);
                    S(:,:,i2) = double(l<p);
                    S(:,:,i2) = ifftshift(S(:,:,i2));
                end
            end
            
            if trajectory==C.LOW_PASS
                S = zeros(phantom_size);
                p = round((n-round(sqrt(n*n*coverage)))/2);
                for i2=1:n
                    S(p:n-p,p:n-p,i2) = 1;
                    S(:,:,i2) = ifftshift(S(:,:,i2));
                end
            end
            
            if trajectory==C.COMPLETE
                S = ones(phantom_size);
            end
            
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        function [ out ] = dct3( in )
            % assuming input to be a 3-D volume of size NxMxP
            M1 = dctmtx(size(in,1));
            M2 = dctmtx(size(in,2));
            M3 = dctmtx(size(in,3));
            out = zeros(size(in));
            for i=1:size(in,3)
                out(:,:,i) = M1 * in(:,:,i) * M2';
            end
            for i=1:size(in,1)
                for j=1:size(in,2)
                    out(i,j,:) = squeeze(out(i,j,:))' * M3';
                end
            end
        end
        
        function [ out ] = idct3( in )
            % assuming input to be a 3-D volume of size NxMxP
            M1 = dctmtx(size(in,1));
            M2 = dctmtx(size(in,2));
            M3 = dctmtx(size(in,3));
            out = zeros(size(in));
            for i=1:size(in,3)
                out(:,:,i) = M1' * in(:,:,i) * M2;
            end
            for i=1:size(in,1)
                for j=1:size(in,2)
                    out(i,j,:) = squeeze(out(i,j,:))' * M3;
                end
            end
        end
        
        function [ theta ] = msfft2( phantom )
            theta = zeros(size(phantom));
            for i=1:size(phantom,3)
                theta(:,:,i) = fft2(phantom(:,:,i));
            end
        end
        
        function [ phantom ] = imsfft2( theta )
            phantom = zeros(size(theta));
            for i=1:size(theta,3)
                phantom(:,:,i) = ifft2(theta(:,:,i));
            end
        end
        
    end
end

function [ S ] = radial( trajectory, n, scale, i2, rot_deg, line_nbr, line_std )
    S = zeros(n);
    L = round(scale); % number of radial lines in the Fourier domain
    for i3=1:line_nbr
        if trajectory==1
            % aperture encompassing all scanning angles (aperture<pi is limited angle)
            aperture  = (pi/180)*180;
            % direction of the scanning beam (middle axis)
            direction = (pi/180)*0+(i2-1)*rot_deg+(i3-1)*2*pi/line_nbr;
        else
            aperture=(pi/180)*90;
            direction=(pi/180)*45+(i2-1)*rot_deg+(i3-1)*2*pi/line_nbr;
        end
        if (pi-aperture)>(aperture/L)
            thc = linspace(-direction-aperture/2, -direction+aperture/2, L);
        else
            thc = linspace(-direction-pi/2, -direction+pi/2-pi/L, L);
        end
        thc = mod(thc,pi);
        for ll = 1:L
            if ((thc(ll) <= pi/4) || (thc(ll) > 3*pi/4))
                yr = round(tan(thc(ll))*(-n/2+1:n/2-1)+n/2+1);
                yr = round(yr + line_std*randn(size(yr)));
                for nn = 1:n-1
                    yr(nn) = max(1,min(yr(nn),n));
                    S(yr(nn),nn+1) = 1;
                end
            else
                xc = round(cot(thc(ll))*(-n/2+1:n/2-1)+n/2+1);
                xc = round(xc + line_std*randn(size(xc)));
                for nn = 1:n-1
                    xc(nn) = max(1,min(xc(nn),n));
                    S(nn+1,xc(nn)) = 1;
                end
            end
        end
    end
end

function [ S ] = spiral( t, n, scale, i2, rot_deg, line_nbr, line_std )
    S = zeros(n);
    for i3=1:line_nbr
        temp = zeros(n);
        x = scale*t .* cos( t+(i2-1)*rot_deg+(i3-1)*2*pi/line_nbr );
        y = scale*t .* sin( t+(i2-1)*rot_deg+(i3-1)*2*pi/line_nbr );
        for i1=1:length(x)
            if round(x(i1)+n/2)<1 || round(y(i1)+n/2)<1 || ...
                    round(x(i1)+n/2)>n || round(y(i1)+n/2)>n
                break
            end
            xr = max(1,min(n,round(x(i1)+n/2)));
            yr = max(1,min(n,round(y(i1)+n/2)));
            if temp(yr,xr)==0
                temp(yr,xr) = 1;
                xr = max(1,min(n,round(x(i1)+n/2+line_std*randn)));
                yr = max(1,min(n,round(y(i1)+n/2+line_std*randn)));
                S(yr,xr) = 1;
            end
        end
    end
end

function [ S ] = log_spiral( t, n, scale, i2, rot_deg, line_nbr, line_std )
    S = zeros(n);
    for i3=1:line_nbr
        temp = zeros(n);
        x = ( exp(scale*t) .* cos(t+(i2-1)*rot_deg+(i3-1)*2*pi/line_nbr) );
        y = ( exp(scale*t) .* sin(t+(i2-1)*rot_deg+(i3-1)*2*pi/line_nbr) );
        % ensure DC is taken
        S(round(n/2),round(n/2)) = 1;
        temp(round(n/2),round(n/2)) = 1;
        for i1=1:length(x)
            if round(x(i1)+n/2)<1 || round(y(i1)+n/2)<1 || ...
                    round(x(i1)+n/2)>n || round(y(i1)+n/2)>n
                break
            end
            xr = max(1,min(n,round(x(i1)+n/2)));
            yr = max(1,min(n,round(y(i1)+n/2)));
            if temp(yr,xr)==0
                temp(yr,xr) = 1;
                xr = max(1,min(n,round(x(i1)+n/2+line_std*randn)));
                yr = max(1,min(n,round(y(i1)+n/2+line_std*randn)));
                S(yr,xr) = 1;
            end
        end
    end
end

function [C] = self_constantsSparseTraj3D
    C = helper.constantsSparseTraj3D;
end

