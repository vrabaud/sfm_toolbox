function anim = addNoise( anim, varargin )
% add noise to the animation
%
% USAGE
%  anim= anim.addNoise()
%
% INPUTS
%  anim     - Animation object (help Animation for details)
%  varargin  - list of paramaters in quotes alternating with their values
%        'noiseS'  - standard variation of the Gaussian noise to add,
%                     if in quotes (e.g. '20'), it is the percentage as
%                     defined in Torresani PAMI08 p. 884
%        'noiseW'  - standard variation of the Gaussian noise to add
%        'doFillS'  - update S with its noisy value: in case of a rigid
%                     scene it will save all the noisy S in a
%                     [ 3 x nPoint x nFrame ]
%        'doFillW'  - re-compute the reprojection values
%
% OUTPUTS
%  anim     - modified Animation
%
% EXAMPLE
%
% See also ANIMATION
%
% Vincent's Structure From Motion Toolbox      Version 3.1
% Copyright (C) 2008-2011 Vincent Rabaud.  [vrabaud-at-cs.ucsd.edu]
% Please email me if you find bugs, or have suggestions or questions!
% Licensed under the GPL [see external/gpl.txt]

[ noiseS noiseW doFillS doFillW ] = getPrmDflt( varargin, {'noiseS',0,...
  'noiseW', 0, 'doFillS',false, 'doFillW',true}, 1);

% deal with the noise in S
if ischar(noiseS)
  if isempty(anim.mask)
    avNormW=norm(anim.W(:),'fro')/anim.nPoint/anim.nFrame;
  else
    avNormW=norm(anim.W(~isnan(anim.W)),'fro')/nnz(anim.mask);
  end
  noiseS=sqrt(str2double(noiseS)/100*avNormW);
end
if noiseS > 0
  if size(anim.S,3)==1
    S = bsxfun(@plus, anim.S, randn( 3, anim.nPoint, anim.nFrame)*noiseS );
  else
    S = anim.S + randn( 3, anim.nPoint, anim.nFrame )*noiseS;
  end
  
  if doFillS; anim.S = S; end
  if doFillW
    SOri = anim.S; anim.S = S; anim.W=generateW(anim); anim.S=SOri;
  end
end

% deal with the noise in W
if noiseW > 0
  W = anim.W + randn( 2, anim.nPoint, anim.nFrame )*noiseW;
  if doFillW, anim.W = W; end
end
