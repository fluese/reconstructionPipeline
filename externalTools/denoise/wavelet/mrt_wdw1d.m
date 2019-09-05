%
% +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% MRT_MATPV 
% +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
%
% MRT_WDW1D MRI processing window
% 
% mrt_wdw1d(name,pts,width,shifts)
% returns a window function of size (pts,1).
%
% input: 
% name      ::  'gauss' 'tukey' 'no'
% pts       ::  number of points 
% width     ::  width, note different definitions of width!!!
%               gauss:  width is defined as the reciprocal of the standard
%                   deviation and is a measure of the width of its FT.
%               tukey:  The 'width' parameter specifies the ratio of taper 
%                   to constant sections. This ratio is normalized to 1 
%                   (i.e., 0 < R < 1).
% shifts    ::  shift relative to center. If omitted, shifts=0. The value 
%               is rounded.
%+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% output:
%           :: vector of size (pts,1) 
%+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
%
% EXAMPLE1:
%           name='gauss'; pts=64; width=4;shifts=-10;
%           w=mrt_wdw1d(name,pts,width,shifts);
%           plot(w)
% EXAMPLE2: 
%           w=mrt_wdw1d('tukey',128,0.7);
%           plot(w)
%
% See also mrt_MATPV, gausswin, tukeywin


%   Reference:
%               no ref.

%   Author(s):  A. Pampel
%               H. Marschner
%   Versions:
%   2012/01/10  AP 

function w=mrt_wdw1d(name,pts,width,shifts)
if(nargin<4);shifts=0;end

shift=abs(round(shifts));

switch lower(name)
    case{'gauss'}
        wt=mygauss(pts+2*shift,width);
        shiftw;
    case{'tukey'}
        wt=mytukey(pts+2*shift,width);
        w=wt(1:pts);
    case{'no'}
        w=ones(pts,1);
    otherwise
        w=zeros(pts,1);
        
end
    function shiftw
        if(sign(shifts)==1)
            w=wt(1:pts);
        else
            w=flipud(wt(1:pts));
        end
    end

    function mg=mygauss(L,a)
        
        error(nargchk(1,2,nargin,'struct'));
        
        % Default value for Alpha
        if nargin < 2 || isempty(a),
            a = 2.5;
        end
        
        % Compute window according to [1]
        N = L-1;
        n = (0:N)'-N/2;
        mg = exp(-(1/2)*(a*n/(N/2)).^2);
    end

    function mt=mytukey(n,r)
        
        error(nargchk(1,2,nargin,'struct'));
        
        % Default value for R parameter.
        if nargin < 2 || isempty(r),
            r = 0.500;
        end
        
        if r <= 0,
            mt = ones(n,1);
        elseif r >= 1,
            mt = hann(n);
        else
            t = linspace(0,1,n)';
            % Defines period of the taper as 1/2 period of a sine wave.
            per = r/2;
            tl = floor(per*(n-1))+1;
            th = n-tl+1;
            % Window is defined in three sections: taper, constant, taper
            mt = [ ((1+cos(pi/per*(t(1:tl) - per)))/2);  ones(th-tl-1,1); ((1+cos(pi/per*(t(th:end) - 1 + per)))/2)];
        end
    end
end