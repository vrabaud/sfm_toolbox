function [l,u,vBest,currBest,lowerBound,isRefined]=bnbRefine(...
  criterionHandle,l,u,vBest,currBest,lowerBound,isRefined)
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
%  criterionHandle - function handle to the criterion
%
% OUTPUTS
%  l             - [ dim x nInterval ] updated l
%  u             - [ dim x nInterval ] updated u
%  vBest         - [ dim x nInterval ] updated vBest
%  currBest      - [ 1 x nInterval ] updated currBest
%  lowerBound    - [ 1 x nInterval ] updated lowerBound
%  isRefined     - [ 1 x nInterval ] updated isRefined
%
% EXAMPLE
%
% See also
%
% Vincent's Structure From Motion Toolbox      Version NEW
% Copyright (C) 2008-2010 Vincent Rabaud.  [vrabaud-at-cs.ucsd.edu]
% Please email me if you find bugs, or have suggestions or questions!
% Licensed under the GPL [see external/gpl.txt]

% figure out the intervals that have not been refined whose lower bound is
% lower than the lowest lowerd bound of the refined intervals
if nnz(isRefined)==0
  [disc,badInd]=min(lowerBound);
else
  badInd=find(lowerBound<min(lowerBound(isRefined)) & ~isRefined);
end

% add the minimum in case it was not refined
[disc,ind]=min(lowerBound);
if ~isRefined(ind); badInd=unique([badInd,ind]); end
minCurrBest=min(currBest);
for ind=badInd
  % make sure the best value is in the interval
  vBestTmp=min([max([vBest(:,ind),l(:,ind)],[],2),u(:,ind)],[],2);

  % check for a few random values in the interval and check if they are
  % better than the value at the optimal v
  allV=[ vBestTmp, bsxfun(@plus,bsxfun(@times,rand(length(l(:,ind)),...
    10000),u(:,ind)-l(:,ind)),l(:,ind))];
  allMin=criterionHandle(allV);
  [ disc, bestInd ]=min(allMin);
  vBestTmp=allV(:,bestInd);

  % perform gradient descent from the best v to find a better solution
  vBestTmp=fmincon(criterionHandle,vBestTmp,[],[],[],[],...
    l(:,ind),u(:,ind),[],optimset(...
    'GradObj','off','Hessian','off','Algorithm','active-set',...
    'Display','off'));
  % make sure the best value is in the interval
  vBestTmp=min([max([vBestTmp,l(:,ind)],[],2),u(:,ind)],[],2);
  % keep that v if the criterion is better than the current one
  currBestTmp=criterionHandle(vBestTmp);

  if currBestTmp<currBest(ind)
    currBest(ind)=currBestTmp; vBest(:,ind)=vBestTmp;
  end
  isRefined(ind)=true;
end

% remove intervals for which the lower bound is higher than the current
% best of another interval
badInterval=find(lowerBound>min(currBest));

if ~isempty(badInterval) && length(badInterval)~=length(currBest)
  l(:,badInterval)=[]; u(:,badInterval)=[]; vBest(:,badInterval)=[];
  lowerBound(badInterval)=[]; currBest(badInterval)=[];
  isRefined(badInterval)=[];
end
