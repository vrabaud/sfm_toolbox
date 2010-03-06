function anim=generateCamFromRt( anim )
% Generate camera position from R and t in an Animation
%
% USAGE
%  anim = anim.generateCamFromRt()
%
% INPUTS
%  anim     - Animation object (help Animation for details)
%
% OUTPUTS
%  anim     - modified Animation object
%
% EXAMPLE
%
% See also GENERATETOYANIMATION
%
% Vincent's Structure From Motion Toolbox      Version NEW
% Copyright (C) 2009 Vincent Rabaud.  [vrabaud-at-cs.ucsd.edu]
% Please email me if you find bugs, or have suggestions or questions!
% Licensed under the Lesser GPL [see external/lgpl.txt]

anim.cam=multiTimes(-anim.R,anim.t,3.1);
