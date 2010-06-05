function [l,u,vBest,currBest,lowerBound,newInd]=bnbBranch(l,u,vBest,...
  currBest,lowerBound)
% Perform the simple branching in brach an d bound
%
% Used in affine and metric Upgrade
%
% USAGE
%   [ l, ul, vBest, currBest, lowerBound, goodInd ]=bnbBranch( l, ...
%    u, vBest, currBest, lowerBound)
%
% INPUTS
%  l             - [ dim x nInterval ] lower bounds
%  u             - [ dim x nInterval ] upper bounds
%  vBest         - [ dim x nInterval ] best value of the parameter in the
%                                      interval
%  currBest      - [ 1 x nInterval ] best value of the criterion on the
%                                      interval
%  lowerBound    - [ 1 x nInterval ] underestimator of the criterion on the
%                                      interval
%
% OUTPUTS
%  l             - [ dim x nInterval ] updated l
%  u             - [ dim x nInterval ] updated u
%  vBest         - [ dim x nInterval ] updated vBest
%  currBest      - [ 1 x nInterval ] updated currBest
%  lowerBound    - [ 1 x nInterval ] updated lowerBound
%  newInd        - array of the new intervals in l,u,vBest ...
%
% EXAMPLE
%
% See also
%
% Vincent's Structure From Motion Toolbox      Version NEW
% Copyright (C) 2008-2010 Vincent Rabaud.  [vrabaud-at-cs.ucsd.edu]
% Please email me if you find bugs, or have suggestions or questions!
% Licensed under the GPL [see external/gpl.txt]

% split the rectangle with the lowest lower bound and the one with the
% lowest current best: both have the best chances to be the right
% interval
[ disc, ind1 ]=min(lowerBound);
[ disc, ind2 ]=min(currBest);
% only do the one with the current best once in a while: if it is indeed
% the best, we are going to refine that interval for no reason
if rand()>0.9; goodInd=unique([ind1,ind2]);
else; goodInd=unique([ind1,ind1]); end

for ind=goodInd(end:-1:1)
  % define the boundaries of the new interval
  [ disc, dim ]=max(u(:,ind)-l(:,ind));
  l(:,end+1)=l(:,ind);
  u(:,end+1)=u(:,ind); u(dim,end)=(u(dim,ind)+l(dim,ind))/2;
  l(:,end+1)=l(:,ind); l(dim,end)=(u(dim,ind)+l(dim,ind))/2;
  u(:,end+1)=u(:,ind);
  vBest(:,end+1:end+2)=[vBest(:,ind),vBest(:,ind)];

  % keep track of the current best and lower bound
  % those can be inherited for the interval we are splitting
  if vBest(dim,ind)>=l(dim,end)
    currBest(end+1:end+2)=[Inf,currBest(ind)];
    lowerBound(end+1:end+2)=[Inf,lowerBound(ind)];
  else
	currBest(end+1:end+2)=[currBest(ind),Inf];
    lowerBound(end+1:end+2)=[lowerBound(ind),Inf];
  end
  % remove the interval we have split
  l(:,ind)=[]; u(:,ind)=[]; vBest(:,ind)=[]; currBest(ind)=[];
  lowerBound(ind)=[];
end

newInd=[max([1,size(l,2)-2*length(goodInd)+1]):size(l,2)];
