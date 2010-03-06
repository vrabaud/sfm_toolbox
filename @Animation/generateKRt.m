function anim = generateKRt(anim,isKIdentity)
% Compute K, R and t using the projection matrices
%
% USAGE
%  anim = anim.generateKRt()
%
% INPUTS
%  anim        - Animation object
%  isKIdentity - flag to tell if is identity or not
%
% OUTPUTS
%  anim       - modified Animation object
%
% EXAMPLE
%
% See also
%
% Vincent's Structure From Motion Toolbox      Version NEW
% Copyright (C) 2009 Vincent Rabaud.  [vrabaud-at-cs.ucsd.edu]
% Please email me if you find bugs, or have suggestions or questions!
% Licensed under the Lesser GPL [see external/lgpl.txt]

if nargin==1 isKIdentity=false; end

anim.R = zeros(3,3,anim.nFrame); anim.t = zeros(3,anim.nFrame);

% fill the different components
if isKIdentity
  for i = 1 : anim.nFrame
    [ anim.R(:,:,i) anim.t(:,i) ] = extractFromP(anim.P(:,:,i),...
      anim.isProj, true);
  end
else
  if anim.isProj anim.K = zeros(5,anim.nFrame);
  else anim.K = zeros(3,anim.nFrame);
  end
  for i = 1 : anim.nFrame
    [ K anim.R(:,:,i) anim.t(:,i) ] = extractFromP(anim.P(:,:,i),
      anim.isProj, false);
    if anim.isProj anim.K(:,i) = [ K(1,1); K(1:2,2); K(1:2,3) ];
    else  anim.K(:,i) = [ K(1,1); K(1:2,2) ]; end
  end
end
