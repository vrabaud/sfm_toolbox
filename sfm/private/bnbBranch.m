function [l,u,vBest,currBest,lowerBound,isRefined,newInd]=bnbBranch(l,u,...
  vBest,currBest,lowerBound,isRefined)
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
%  isRefined     - [ 1 x nInterval ] if true, the result in there has been
%                                      refined
%
% OUTPUTS
%  l             - [ dim x nInterval ] updated l
%  u             - [ dim x nInterval ] updated u
%  vBest         - [ dim x nInterval ] updated vBest
%  currBest      - [ 1 x nInterval ] updated currBest
%  lowerBound    - [ 1 x nInterval ] updated lowerBound
%  isRefined     - [ 1 x nInterval ] updated isRefined
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
if rand()>0.9; goodInd=unique([ind1,ind1]);
else; goodInd=unique([ind1,ind1]); end

% also add the biggest interval
[ disc, ind3 ]=max(prod(u-l,1));
%  goodInd(end+1)=ind3;
%  goodInd(end+1)=1;

goodInd=unique(goodInd);

for ind=goodInd(end:-1:1)
  % define the boundaries of the new interval
  % one way of doing it is by simply splitting in the middle of the largest
  % dimension


  % you can also split where there is already a minimum, if it is not too
  % far from the center
  [ disc, dim ]=max(u(:,ind)-l(:,ind));
  if 0 && abs(vBest(dim,ind)-(u(dim,ind)+l(dim,ind))/2) < 0.4*(u(dim,ind)-l(dim,ind))/2
    l(:,end+1)=l(:,ind);
    u(:,end+1)=u(:,ind); u(dim,end)=vBest(dim,ind);
    l(:,end+1)=l(:,ind); l(dim,end)=vBest(dim,ind);
    u(:,end+1)=u(:,ind);
  else
    [ disc, dim ]=max(u(:,ind)-l(:,ind));
    l=l(:,[1:end,ind,ind]); u=u(:,[1:end,ind,ind]);
    u(dim,end-1)=(u(dim,ind)+l(dim,ind))/2;
    l(dim,end)=(u(dim,ind)+l(dim,ind))/2;
  end

  vBest(:,end+1:end+2)=[vBest(:,ind),vBest(:,ind)];

  % keep track of the current best and lower bound
  % those can be inherited from the interval we are splitting
  if vBest(dim,ind)>=l(dim,end)
    currBest(end+1:end+2)=[Inf,currBest(ind)];
    lowerBound(end+1:end+2)=[Inf,lowerBound(ind)];
    isRefined(end+1:end+2)=[false,isRefined(ind)];
  else
    currBest(end+1:end+2)=[currBest(ind),Inf];
    lowerBound(end+1:end+2)=[lowerBound(ind),Inf];
    isRefined(end+1:end+2)=[isRefined(ind),false];
  end
  % remove the interval we have split
  l(:,ind)=[]; u(:,ind)=[]; vBest(:,ind)=[]; currBest(ind)=[];
  lowerBound(ind)=[]; isRefined(ind)=[];
end

newInd=[size(l,2)-2*length(goodInd)+1:size(l,2)];
