function AverageFactors=getAverageFactors(twix_obj, Lines, Partitions)
% Function to get the weighting factors
% 
% ************************************************************************
% Version 1.0                                                  15.02.2019
% Falk Luesebrink              falk dot luesebrink at med dot ovgu dot de
% ************************************************************************

% get the data points
P = [twix_obj.image.Lin; twix_obj.image.Par];

% initialize map
a = zeros(twix_obj.image.NLin, twix_obj.image.NPar);

% get average factors
A = twix_obj.image.freeParam(1,:);

% assign weighting factors to map
for j = 1:size(P,2)
    a(P(1,j),P(2,j)) = A(j);
end

% correct fo zeros
m = (a==0);
a(m)=1;

% compute AverageFactors
AverageFactors = a( Lines(1):Lines(end), Partitions(1):Partitions(end) );
