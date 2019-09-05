function [wdwnd,parchecked]=mrt_wdw(par, sz)
%
% +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% MRT_MATPV
% +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
%
% [wdwnd,parchecked]=mrt_wdwnd(para)
% returns a 1D,2D, or 3D MRI window
%
% input:
%
% a) no input
%       for testing, a standard para set is set
%
% b) structures containing the fields 'name', 'size' and 'width', 'shifts'
%       'shifts' can be omitted and will be set to zero.
%
%    See mrt_wdw1d for definitions
%
% output:
%    [wdwnd,parchecked]
%           wdwnd :: 1D,2D,3D, matrix
%           parchecked :: the para structure that was used for creating
%           wdwnd
% EXAMPLE1:
%           para.name={'tukey','gauss','gauss'};
%           para.size=[128 64 256];
%           para.width=[0.5 2 2];
%           para.shifts=[0 2 -3];
%           [w,pp]=mrt_wdw(para);
%           imagesc(squeeze(w(:,:,64)))
% EXAMPLE2:
%           para.name={'tukey'};
%           para.size=[128];
%           para.width=[0.5];
%           para.shifts=[0];
%           [w,pp]=mrt_wdw(para);
%           plot(w)
%
%
% See also mrt_MATPV, mrt_wdw1d, gausswin, tukeywin


%   Reference:
%               no ref.

%   Author(s):  A. Pampel
%               H. Marschner
%   Versions:
%   2012/01/10  AP

% if (nargin<1 || isstruct(par)==0)
%     % For testing only!
%     % Wrong in reco !!!!!!
%     par.name={'tukey' 'gauss' 'tukey'};
%     par.size=[128 256 64];
%     par.width=[0.5 2.5 0.5];
%     par.shifts=[0 0 0];
% elseif (nargin==1 && isstruct(par)==1)
%     repaired=check_par;
%     if(~repaired);
%         disp('window set to one')
%         wdwnd=1;
%         parchecked='ERROOOOOOOR';
%         return;
%     end
% else
%     error('unknown error! See help.')
% end

name=par.name(:);
si=sz;
wi=par.width(:);
shi=par.shifts(:);

% Reading dimensions.
dname=size(name,1);
dsi=size(si,1);
dwi=size(wi,1);
dshi=size(shi,1);

switch dsi
    case 3
        % check consistency
        if(dname+dsi+dwi+dshi<12)
            error('Window parameter dimensions do not agree!')
        end
        % end of checks
        
        % In case of name only containing 'no'.
        if sum(strcmp({name{1},name{2},name{3}},{'no','no','no'}))==3
            wdwnd=1;
            disp('Window set to 1! No 3D window created!');
        else
            wdw3d=zeros(si');% creation of arrays of final dimensions
            wdw3dt=zeros(si');
            
            w1=mrt_wdw1d(name{1},si(1),wi(1),shi(1));
            w2=mrt_wdw1d(name{2},si(2),wi(2),shi(2));
            w3=mrt_wdw1d(name{3},si(3),wi(3),shi(3));
            wdw2d=w1*w2';
            
            for al=1:si(3)
                wdw3dt(:,:,al)=wdw2d;
            end
            
            w3t=repmat(w3,1,si(1))';
            for ka=1:si(2)
                wdw3d(:,ka,:)=squeeze(wdw3dt(:,ka,:)).*w3t;
            end
            
            wdwnd=wdw3d;
            
        end
        
        disp('3D window');
        
    case 2
        % checks
        if(dname+dsi+dwi+dshi<8)
            error('Wrong dimensions of input data para.wpar! See help.')
        end
        % end of checks
        
        % In case of name only containing 'no'.
        if sum(strcmp({name{1},name{2}},{'no','no'}))==2
            wdwnd=1;
        else
            w1=mrt_wdw1d(name{1},si(1),wi(1),shi(1));
            w2=mrt_wdw1d(name{2},si(2),wi(2),shi(2));
            
            wdwnd=w1*w2';
        end
        
        disp('2D window');
        
    case 1
        if(strcmpi(name,{'NO'})) 
            wdwnd=1;
        else
            wdwnd=mrt_wdw1d(name{1},si(1),wi(1),shi(1));
        end
    otherwise
        error('These are no 1, 2 or 3 dimensional data!')
end

parchecked=par;

%+++++++++++++++++++++++++++++++++++++++++++++++++++++
% NESTED FUNCTIONS
%+++++++++++++++++++++++++++++++++++++++++++++++++++++

    function repaired=check_par
        checkforthisfields={'name','size','width'};
        field=isfield(par,checkforthisfields);
        
        
        if(sum(field)~=3)
            repaired=false;
            disp(['Inputerror! --> ' ...
                'Field(s) missing : ' ....
                cell2mat(checkforthisfields(~field)) ])
            return
        end
        
        if(~isequal(length(par.size(:)),length(par.width(:)),length(par.name(:))))
            repaired=false;
            disp('Window parameter dimensions do not agree!')
            return
        end
        
        if(isfield(par,{'shifts'}) && isequal(length(par.size(:)),length(par.shifts(:))))
            disp('window par seems ok...')
            repaired=true;
        elseif(~isfield(par,{'shifts'}))
            par.shifts=zeros(size(par.size));
            disp('window par seems ok...')
            disp('all shifts set to zero !!!')
            repaired=true;
        else
            repaired=false;
            disp('Inputerror! Cannot repair input.')
        end
    end
end
