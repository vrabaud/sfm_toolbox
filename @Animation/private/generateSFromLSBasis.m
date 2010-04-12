function anim = generateSFromLSBasis( anim )
% Generate 3D data from the basis of an Animation
%
% USAGE
%  anim = anim.generateSFromLSBasis()
%
% INPUTS
%  anim     - Animation object (help Animation for details)
%
% OUTPUTS
%  anim     - modified Animation anime
%
% EXAMPLE
%
% See also GENERATETOYANIMATION
%
% Vincent's Structure From Motion Toolbox      Version 3.0
% Copyright (C) 2009 Vincent Rabaud.  [vrabaud-at-cs.ucsd.edu]
% Please email me if you find bugs, or have suggestions or questions!
% Licensed under the GPL [see external/gpl.txt]

if size(anim.l,1)>=(size(anim.SBasis,3)-1) && ~isempty(anim.SBasis)
  anim.nPoint=size(anim.SBasis,2);
  switch size(anim.l,1)
    case anim.nBasis,
      anim.S=reshape(sum(bsxfun(@times,reshape(anim.l,1,1,anim.nBasis,...
        anim.nFrame),anim.SBasis),3),3,anim.nPoint,anim.nFrame);
    case anim.nBasis-1,
      % add the first basis with a coeff of 1
      anim.S=bsxfun(@plus,anim.SBasis(:,:,1), reshape(sum(bsxfun(@times,...
        reshape(anim.l,1,1,anim.nBasis-1,anim.nFrame),...
        anim.SBasis(:,:,2:end)),3),3,anim.nPoint,anim.nFrame));
    otherwise,
      error('l and nBasis do not have matching dimensions');
  end
end
