function S=generateSAbsolute( anim )
% Generate S tranformed by R and t
%
% USAGE
%  anim = anim.generateSAbsolute()
%
% INPUTS
%  anim     - Animation object (help Animation for details)
%
% OUTPUTS
%  S        - anim.S after R, t have been applied (or homographies)
%
% EXAMPLE
%
% See also GENERATETOYANIMATION
%
% Vincent's Structure From Motion Toolbox      Version 3.0
% Copyright (C) 2008-2011 Vincent Rabaud.  [vrabaud-at-cs.ucsd.edu]
% Please email me if you find bugs, or have suggestions or questions!
% Licensed under the GPL [see external/gpl.txt]

if isempty(anim.S)
  warning('S is not defined nor computable'); S=[]; return;
end
if (isempty(anim.R) || isempty(anim.t)) && (isempty(anim.P))
  warning('Need to define R,t or P'); S=[]; return;
end


if size(anim.P,2)==3
  % homography case
  if size(anim.S,3)==1
    S = multiTimes( anim.P, anim.S, 1 );
  else
    S = multiTimes( anim.P, anim.S, 2 );
  end
else
  % compute the projected points for normal camera and 3D points
  if isempty(anim.R)
    if size(anim.S,3)==1
      S = bsxfun( @plus, multiTimes( anim.P(:,1:3,:), anim.S, 1 ), ...
        anim.P(:,4,:) );
    else
      S = bsxfun( @plus, multiTimes( anim.P(:,1:3,:), anim.S, 2 ), ...
        anim.P(:,4,:) );
    end
  else
    if size(anim.S,3)==1
      S = bsxfun( @plus, multiTimes( anim.R, anim.S, 1 ), ...
        reshape(anim.t, 3, 1, anim.nFrame ) );
    else
      S = bsxfun( @plus, multiTimes( anim.R, anim.S, 2 ), ...
        reshape(anim.t, 3, 1, anim.nFrame ) );
    end
  end
end
