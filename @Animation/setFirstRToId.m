function [anim,HEye] = setFirstPRtToId( anim )
% Sets the first P or (R,t) of an Animation to eye(3,4)
%
% The function changes all the S and SBasis accordingly to preserve the
% projections.
%
% if R and t is defined, the rigid transform H is computed such that
%   [R(:,:,1),t(:,1)]*H=eye(3,4)
%   H is of the form [rotation,translation; 0 0 0 1]
%
% if R is not defined but P is, the homography H is computed such that:
%   P(:,:,1)*H=eye(3,4)
%
% USAGE
%  anim = anim.setFirstRToId()
%
% INPUTS
%  anim     - Animation object (help Animation for details)
%
% OUTPUTS
%  anim     - modified Animation
%  HEye     - [3x4] rigid transform if (R,t) are defined, or [3x4]
%             projective transform is P is defined
%
% EXAMPLE
%
% See also GENERATETOYANIMATION
%
% Vincent's Structure From Motion Toolbox      Version 3.0
% Copyright (C) 2009 Vincent Rabaud.  [vrabaud-at-cs.ucsd.edu]
% Please email me if you find bugs, or have suggestions or questions!
% Licensed under the GPL [see external/gpl.txt]

if isempty(anim.R)
  HEye=anim.P(:,:,1)\eye(3,3);
  [disc,disc,V]=svd(anim.P(:,:,1));
  HEye(:,4)=V(:,4);
  % apply the homography to P
  anim.P=multiTimes(anim.P,HEye,1);
else
  HEye=anim.R(:,:,1)';
  HEye(:,4)=-anim.R(:,:,1)'*anim.t(:,1);
  % Re-generate the rotations
  % need the following line for compatibility with octave ... otherwise
  % I would use end
  PAll=anim.P;
  if anim.nFrame>=2
    PAll(:,:,2:anim.nFrame)=multiTimes(PAll(:,:,2:anim.nFrame),HEye,1);
  end
  PAll(:,:,1)=eye(3,4); % Just for numerical stability

  % assign values to R and t
  anim.R=anim.P(:,1:3,:); anim.t=reshape(anim.P(:,4,:),3,anim.nFrame);
end

if anim.nBasis~=0
  % Re-generate the basis
  % use subsasgn so that S is modified at the same time too
  SBasis=normalizePoint(multiTimes(inv(HEye),normalizePoint(anim.SBasis,-4),1.2),4);
  anim=subsasgn(anim,struct('type','.','subs','SBasis'),SBasis);
else
  anim.S=normalizePoint(multiTimes(inv(HEye),normalizePoint(anim.S,-4),1.2),4);
end
