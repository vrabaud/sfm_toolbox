function [ grSOrtho, H ] = msfmGradientLsml( anim, d )
% Compute the LSML gradient in MSFM from Rabaud CVPR08
%
% Just check Rabaud's CVPR08 paper
% the notations follow Rabaud's PhD thesis: Appendix D3
%
% USAGE
%  [ grSOrtho, H ] = msfmGradientLsml( anim, d )
%
% INPUTS
%  anim          - Animation object
%  d             - dimensionality of the shape manifold
%
% OUTPUTS
%  grSOrtho      - LSML gradient, orthogonal to the manifold
%  H             - [ 3nPoint x d x nFrame ] tangent basis at every point
%                  of the manifold
%
% EXAMPLE
%
% See also MSFMGRADIENTSRT
%
% Vincent's Structure From Motion Toolbox      Version 3.0
% Copyright (C) 2008-2011 Vincent Rabaud.  [vrabaud-at-cs.ucsd.edu]
% Please email me if you find bugs, or have suggestions or questions!
% Licensed under the GPL [see external/gpl.txt]

nFrame = anim.nFrame; nPoint = anim.nPoint;

X3 = reshape(anim.S,3*nPoint, nFrame);

% figure out neighbors in the data we have
neigh = computeNeighbor( X3,'k',3,'forceConn',2);
Nmat = neigh.Nmat;

% force temporally closeby frames to be considered together
row = zeros(1,nFrame); row(2) = 1;
Nmat(logical(toeplitz(row))) = 1;

% learn the manifold from here
pTh=struct('d',d,'rbfK',5,'nRestart',5,'nSamples',250,'show',0);
manifold=lsmlOptimizeTh(X3,Nmat,pTh);

% get the gradient for denoising
[Chi3, Chi] = lsmlDenoise( X3, Nmat, manifold, 'nItr', 1 );

grSOrtho = reshape(X3-Chi,3,nPoint,nFrame);

if nargout >=2; H = lsmlComputeH( X3, manifold, 1 ); end
end
