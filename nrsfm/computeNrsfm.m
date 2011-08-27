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
%       - 'nItrSBA'       - number of bundle adjustment iterations
%     for THB, PAMI08
%       - 'nItr'          - maximum number of EM iterations (nItr)
%     for MSFM
%     for CSFM
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
% Vincent's Structure From Motion Toolbox      Version NEW
% Copyright (C) 2008-2011 Vincent Rabaud.  [vrabaud-at-cs.ucsd.edu]
% Please email me if you find bugs, or have suggestions or questions!
% Licensed under the GPL [see external/gpl.txt]

[ nBasis nItr nTriplet nPair nItrSBA ] = ...
  getPrmDflt( varargin, { 'nBasis' [] 'nItr' 50 'nTriplet' 0 'nPair' 0 ...
  'nItrSBA' 0}, 1 );

nFrame = size( W, 3 ); nPoint = size( W, 2 );
anim=Animation('W',W);
hasAnyNan=any(isnan(W));

switch method
  case 1
    % 2D motion resulting from orthographic projection (Eq (1))
    p2_obs = reshape( permute( W, [ 3 1 2 ] ), [], nPoint );
    
    % runs the non-rigid structure from motion algorithm
    use_lds = 0;
    tol = 0.0001;
    
    % build the mask of themissing data
    if ~isempty(anim.mask); MD = ~anim.mask;
    else MD = zeros( nFrame, nPoint );
    end
    
    [ disc SBar V R t Z ] = em_sfm(p2_obs, MD, nBasis-1, use_lds, ...
      tol, nItr);
    
    anim.isProj = false;
    anim.W = W;
    RR=zeros(3,3,nFrame); tt=zeros(3,nFrame);
    for i=1:nFrame; RR(:,:,i) = R{i}; tt(:,i) = R{i}*[t(i,:)';0]; end
    anim.R=RR; anim.t=tt;
    if nBasis==1; anim.S=SBar; anim.l=[]; anim.SBasis=[];
    else
      anim.SBasis=SBar;
      for i=1:nBasis-1; anim.SBasis(:,:,i+1)=V(3*i-2:3*i,:); end
      anim.l=Z';
    end
  case 2 % Xiao Kanade
    if hasAnyNan
      warning('Method not supported with missing entries');
      return;
    end
    if exist('OCTAVE_VERSION','builtin')==5
      warning('Method not supported under Octave');
      return;
    end
    anim = nrsfmXiaoKanade( W, nBasis );
  case 3,
    if hasAnyNan
      warning('Method not supported with missing entries');
      return;
    end
    if exist('OCTAVE_VERSION','builtin')==5
      warning('Method not supported under Octave');
      return;
    end
    anim = nrsfmMsfm( W, nBasis, nTriplet, nPair, nItr );
  case 4,
    if hasAnyNan
      warning('Method not supported with missing entries');
      return;
    end
    if exist('OCTAVE_VERSION','builtin')==5
      warning('Method not supported under Octave');
      return;
    end
    anim = nrsfmCsfm( W, nBasis, nTriplet, nPair, nItr );
end
% Perform gradient descent on the different coefficients
if nItrSBA>0 && (method==1 || method==2)
  anim = bundleAdjustment(anim,'nItr',nItrSBA);
end
