function theta=stein(y,sigmanoise,b)
% James–Stein estimator
% https://en.wikipedia.org/wiki/James%E2%80%93Stein_estimator

n=numel(y);
if(nargin<3)
    b=zeros(size(y));
end

% Formula from Kockelkorn "Statistik" identical ????
%mu=(1-((n-2)*sigmanoise^2)/(norm(y-b,2).^2))*y+(((n-2)*sigmanoise^2)/(norm(y-b,2).^2))*b;
%theta=mu 

klamm=1-((n-3)*sigmanoise^2)/(norm(y-b,2).^2);
klamm(klamm<0)=0;
theta=klamm*(y-b)+b;
