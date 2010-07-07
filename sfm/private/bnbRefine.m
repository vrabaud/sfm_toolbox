function [vBest,currBest]=bnbRefine(l,u,vBest,criterionHandle)
% Refine a result after bounding has been performed
%
% Used in affine and metric Upgrade
%
% USAGE
%   [vBest,currBest]=bnbBranch(l,u,vBest,criterionHandle)
%
% INPUTS
%  l               - [ dim x 1 ] lower bounds
%  u               - [ dim x 1 ] upper bounds
%  vBest           - [ dim x 1 ] best value of the parameter in the
%                                      interval
%  criterionHandle - function handle to the criterion
%
% OUTPUTS
%  vBest         - [ dim x 1 ] refined vBest
%  currBest      - best criterion value
%
% EXAMPLE
%
% See also
%
% Vincent's Structure From Motion Toolbox      Version NEW
% Copyright (C) 2008-2010 Vincent Rabaud.  [vrabaud-at-cs.ucsd.edu]
% Please email me if         you find bugs, or have suggestions or questions!
% Licensed under the GPL [see external/gpl.txt]

% make sure the best value is in the interval
vBest=min([max([vBest,l],[],2),u],[],2);

% check for a few random values in the interval and check if they are
% better than the value at the optimal v
allV=[ vBest, bsxfun(@plus,bsxfun(@times,rand(length(l),10000),u-l),l)];
allMin=criterionHandle(allV);
[ disc, bestInd ]=min(allMin);
vBest=allV(:,bestInd);

% perform gradient descent from the best v to find a better solution
vBest=fmincon(criterionHandle,vBest,[],[],[],[],l,u,[],optimset(...
  'GradObj','off','Hessian','off','Algorithm','active-set',...
  'Display','off'));
% make sure the best value is in the interval
vBest=min([max([vBest,l],[],2),u],[],2);
% keep that v if the criterion is better than the current one
currBest=criterionHandle(vBest);
