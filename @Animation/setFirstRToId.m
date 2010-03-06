function anim = setFirstRToId( anim )
% Sets the first rotation of an Animation to Id
%
% Changes all the R, S and SBasis accordingly
%
% USAGE
%  anim = anim.setFirstRToId()
%
% INPUTS
%  anim     - Animation object (help Animation for details)
%
% OUTPUTS
%  anim     - modified Animation
%
% EXAMPLE
%
% See also GENERATETOYANIMATION
%
% Vincent's Structure From Motion Toolbox      Version NEW
% Copyright (C) 2009 Vincent Rabaud.  [vrabaud-at-cs.ucsd.edu]
% Please email me if you find bugs, or have suggestions or questions!
% Licensed under the Lesser GPL [see external/lgpl.txt]

R=anim.R(:,:,1);
% Re-generate the rotations
% need the following line forcompatibility with octave ... otherwise
% I would use end
RAll=anim.R;
if anim.nFrame>=2
  RAll(:,:,2:anim.nFrame)=multiTimes( anim.R(:,:,2:anim.nFrame), R', 1 );
end
RAll(:,:,1)=eye(3); % Just for numerical stability

% use subsasgn so that P is modified at the same time too
anim=subsasgn(anim,struct('type','.','subs','R'),RAll);

if anim.nBasis~=0
  % Re-generate the basis
  % use subsasgn so that S is modified at the same time too
  anim=subsasgn(anim,struct('type','.','subs','SBasis'),...
    multiTimes( R, anim.SBasis, 1.2 ));
else
  anim.S = multiTimes( R, anim.S, 1.2 );
end
