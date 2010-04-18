function [ animArray SimMin ] = nrsfmCsfm( W, nBasisArray, nTriplet, nPair, ...
  nItr )
% Compute orthographic non-rigid structure from motion using Rabaud CVPR09
%
% Just check the Rabaud Belongie CVPR09 paper
%
% USAGE
%  [ animArray SimMin ] = nrsfmCsfm( W, nBasisArray, nTriplet, nPair, ...
%  nItr, SimMin )
%
% INPUTS
%  W             - [ 2 x nPoint x nFrame ] set of 2D points
%  nBasisArray   - array of basis numbers to use (useful as some computation
%                  for CSFM is independent from the basis number)
%                  can also be just one number
%  nTriplet      - number of triplets of frames to use for GNMDS constraints
%  nPair         - number of pairs of frames to use for GNMDS constraints
%  nItr          - [50] number of iterations to use for bundle adjustment
%
% OUTPUTS
%  animArray     - array containing one reconstructed anim object per basis
%                  number in nBasisArray. Just an anim if nBasisArray is a
%                  number
%  SimMin        - similarity matrix between the pairs of frames
%
% EXAMPLE
%
% See also
%
% Vincent's Structure From Motion Toolbox      Version 3.0
% Copyright (C) 2009 Vincent Rabaud.  [vrabaud-at-cs.ucsd.edu]
% Please email me if you find bugs, or have suggestions or questions!
% Licensed under the GPL [see external/gpl.txt]

if nargin<=4; nItr=50; end

% pre-processing
anim=Animation('W',W);
[ anim anim.t ]=anim.centerW();
anim.t(3,:)=0;

% Compute similarities within views
SimMin = computeAnimSimilarity( anim, -2, 'doDisplay', false);

% Compute samples to study with GNMDS
minDist = 5;
samp = csfmComputeSample( anim, SimMin, nTriplet, minDist );

% perform gnmds
[ lTot K prob ] = csfmGnmds( samp, SimMin, nTriplet, nPair, ...
  max(nBasisArray), 1e1, 0 );

animArray = cell(1,0);
for nBasis = nBasisArray
  anim.l=lTot(1:nBasis-1,:);
  
  % Compute the best basis knowing the coefficients
  anim.SBasis=csfmComputeSBasis( anim.l, anim.W, 5, 1 );
  % anim.S is updated automatically
  
  % Compute the best rotation matrices
  anim.R = computeOrientation(anim.W,anim.S,'exteriorSequence');
  
  if nItr>0
    anim = bundleAdjustment( anim, 'nItr', nItr );
  end
  
  animArray{end+1} = anim;
  save temp;
end

if length(nBasisArray)==1; animArray = animArray(1); end
