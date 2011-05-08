function [ P, anim ] = generateP(anim, doFillP)
% Compute the projection matrices using K, R and t
%
% USAGE
%  P = anim.generateP()
% [ P anim ] = generateP(anim, true)
%
% INPUTS
%  anim       - Animation object
%  doFillP    - if set to true and anim.P is empty, it will fill anim.P
%
% OUTPUTS 'reproj'
%  P       - [ 3 x 4 x nFrame ] projection matrix
%
% EXAMPLE
%
% See also
%
% Vincent's Structure From Motion Toolbox      Version 3.1
% Copyright (C) 2008-2011 Vincent Rabaud.  [vrabaud-at-cs.ucsd.edu]
% Please email me if you find bugs, or have suggestions or questions!
% Licensed under the GPL [see external/gpl.txt]

if nargin==1 || nargout==1; doFillP = false; end

P = zeros( 3, 4, anim.nFrame );
if anim.isProj
  P(:,1:3,:) = anim.R;
  P(:,4,:) = reshape(anim.t,3,1,anim.nFrame);
else
  P(1:2,1:3,:) = anim.R(1:2,:,:);
  P(1:2,4,:) = reshape(anim.t(1:2,:),2,1,anim.nFrame);
  P(3,4,:) = 1;
end

% Apply calibration matrix if any
if ~isempty(anim.K)
  if anim.isProj
    if size(anim.K,2)==1
      P(1,:,:) = sum(bsxfun(@times,anim.K([1,2,4]),P),1);
      P(2,:,:) = anim.K(3)*P(2,:,:)+anim.K(5)*P(3,:,:);
    else
      P(1,:,:) = multiTimes(anim.K([1,2,4],:),P,3.2);
      P(2,:,:) = multiTimes(anim.K([3,5],:),P(2:3,:,:),3.2);
    end
  else
    if size(anim.K,2)==1
      P(1,:,:) = anim.K(1)*P(1,:,:)+anim.K(2)*P(2,:,:);
      P(2,:,:) = anim.K(3)*P(2,:,:);
    else
      P(1,:,:) = multiTimes(anim.K([1,2],:),P(1:2,:,:),3.2);
      P(2,:,:) = multiTimes(anim.K(3,:),P(2,:,:),3.2);
    end
  end
end

% in the case where we want to save the results
if doFillP
  % use subsasgn so that K,R,t are modified at the same time too
  anim=subsasgn(anim,struct('type','.','subs','P'),P);
end
