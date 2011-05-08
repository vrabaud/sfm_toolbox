function [ anim meanW ]=centerW( anim )
% Generate camera position from R and t in an Animation
%
% USAGE
%  [ anim meanW ]= anim.centerW()
%
% INPUTS
%  anim     - Animation object (help Animation for details)
%
% OUTPUTS
%  anim     - modified Animation object
%  meanW    - mean of the measurements (over the points)
%
% EXAMPLE
%
% See also ANIMATION
%
% Vincent's Structure From Motion Toolbox      Version 3.0
% Copyright (C) 2008-2011 Vincent Rabaud.  [vrabaud-at-cs.ucsd.edu]
% Please email me if you find bugs, or have suggestions or questions!
% Licensed under the GPL [see external/gpl.txt]

meanW=mean(anim.W,2);

anim.W = bsxfun(@minus, anim.W, meanW);

if nargout==2
  meanW=squeeze(meanW);
end
