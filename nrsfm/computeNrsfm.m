function anim = computeNrsfm( method, W, varargin )
% Compute orthographic non-rigid structure from motion
%
% Different methods are implemented:
% - method 1: Torresani (PAMI 08)
% - method 2: Xiao, Kanade (IJCV 06)
% - method 3: MSFM Rabaud, Belongie (CVPR 08)
% - method 4: CSFM Rabaud, Belongie (CVPR 09)
%
% USAGE
%  anim = computeNrsfm( method, W, varargin )
%
% INPUTS
%  W             - [ 2 x nPoint x nFrame ] set of 2D points
%  varargin   - list of paramaters in quotes alternating with their values
%       - 'nBasis'        - number of shape bases.  If Torresani's method,
%                           the first shape has a shape coefficient of 1,
%                           and is therefore considered constant.
%     for THB, PAMI08
%       - 'nItr'          - maximum number of EM iterations (nItr)
%     for MSFM
%     for TSFM
%       - 'nTriplet'      -
%       - 'nPair'
%       - 'nItr'          - number of gradient descent iterations
%
% OUTPUTS
%  anim          - Animation object (help Animation for details)
%
% EXAMPLE
%
% See also COMPUTESMFROMW
%
% Vincent's Structure From Motion Toolbox      Version 2.0\n
% Copyright (C) 2009 Vincent Rabaud.  [vrabaud-at-cs.ucsd.edu]
% Please email me if you find bugs, or have suggestions or questions!
% Licensed under the Lesser GPL [see external/lgpl.txt]

[ nBasis nItr nTriplet nPair ] = ...
  getPrmDflt( varargin, { 'nBasis' [] 'nItr' 50 'nTriplet' 0 'nPair' 0 ...
  }, 1 );

nFrame = size( W, 3 ); nPoint = size( W, 2 );
anim=Animation('W',W);

switch method
  case 1
    % 2D motion resulting from orthographic projection (Eq (1))
    p2_obs = reshape( permute( W, [ 3 1 2 ] ), [], nPoint );

    % runs the non-rigid structure from motion algorithm
    use_lds = 0;

    tol = 0.0001;
    MD = zeros( nFrame, nPoint );

    [ disc SBar V R t Z ] = em_sfm(p2_obs, MD, nBasis-1, use_lds, ...
      tol, nItr);

    anim.isProj = false;
    anim.W = W;
    anim.t(3,:)=0;
    for i=1:nFrame
      anim.R(:,:,i) = R{i};
      anim.t(:,i) = R{i}*[t(i,:)';0];
    end
    if nBasis==1; anim.S=SBar; anim.l=[]; anim.SBasis=[];
    else
      anim.l=Z'; anim.SBasis=SBar;
      for i=1:nBasis-1; anim.SBasis(:,:,i+1)=V(3*i-2:3*i,:); end
    end
  case 2 % Xiao Kanade
    anim = nrsfmXiaoKanade( W, nBasis );
  case 3,
    anim = nrsfmMsfm( W, nBasis, nTriplet, nPair, nItr );
  case 4,
    anim = nrsfmCsfm( W, nBasis, nTriplet, nPair, nItr );
end
% Perform gradient descent on the different coefficients
% [ S P K R t l ] = bundleAdjustment( false, Inf, W, animi.S, anim.R, ...
%   anim.t, 'l', anim.l )
