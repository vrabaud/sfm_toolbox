function anim = computeSMFromW( isProj, W, varargin )
% Recover rigid structure and motion from 2D calibrated measurements W
%
% Essential matrix decomposition
% Reference: HZ2, p259, Result 9.19, and p. 294
% isProj && nFrame==2 && isCalibrated
% fast
%
% Normalized 8-point algorithm
% Reference: HZ2, p279 and alg 11.1 p. 282
% isProj && nFrame==2 && ~isCalibrated
% fast
%
% Projective Factorization
% Reference: HZ2, p445, Algorithm 18.2
% Sturm and Triggs ECCV 96
% A Factorization Based Algorithm for Multi-Image Projective Structure and
% Motion
% isProj && nFrame>2 && method==0
%
% Projective Factorization
% Reference: Iterative Extensions of the Sturm/Triggs Algorithm:
% Convergence and Nonconvergence, from Oliensis, Hartley, PAMI 07
% isProj && nFrame>2 && method==Inf
% slower
%
% Gold Standard for Affine camera matrix
% Reference: HZ2, p351, Algorithm 14.1
% ~isProj && nFrame==2 && onlyErrorFlag
% fast
%
% Tomasi Kanade without the metric constraint
% Affine camera matrix, MLE estimation (Tomasi Kanade)
% Reference: HZ2, p437, Algorithm 18.1
% ~isProj && nFrame>2
% fast
%
% Tomasi Kanade with the metric constraint
% ~isProj && nFrame>2 && isCalibrated && method=0
% fast
%
% Tomasi Kanade with the metric constraint (more accurate computation)
% ~isProj && nFrame>2 && isCalibrated && method=inf
% slow
%
% Tomasi Kanade and autocalibration
% ~isProj && nFrame>2 && ~isCalibrated
% slow
%
% if doAffineUpgrade and/or doMetricUpgrade are true, affine and metric
% upgrades from Chandraker IJCV 2009 is applied
% slow
%
% If there are any missing entries:
% Low-Rank Matrix Fitting Based on Subspace Perturbation Analysis
% with Applications to Structure from Motion
% Hongjun Jia, Aleix M. Martinez, PAMI 08
%
% Returns the structure and camera parameters in an Animation object
%
% USAGE
%   anim = computeSMFromW( isProj, W, 'method',method )
%   anim = computeSMFromW( isProj, W,'method',Inf,'onlyErrorFlag',true )
%
% INPUTS
%  isProj     - flag indicating if the camera is projective
%  W          - [] [ 2 x nPoint x nFrame ] 2D projected features (NaN if
%               missing entry
%  varargin   - list of paramaters in quotes alternating with their values
%       - 'isCalibrated' [false] flag indicating if the cameras are
%                        calibrated (intrinsic parameters are identity)
%       - 'doAffineUpgrade' [false] flag indicating if the affine upgrade
%                        is computed
%       - 'doMetricUpgrade' [false] flag indicating if the metric upgrade
%                        is computed
%       - 'K',[] [3 x 1 ], [ 3 x nFrame ] calibration parameters
%                       (or [5 x 1 ], [ 5 x nFrame ] when projective)
%                       If given, they won't be optimized upon, you will
%                       need to run bundle adjustment after that
%       - 'KFull',[] [3 x 3 ] or [ 3 x 3 x nFrame ] calibration matrices
%       - 'method', [inf] method for performing SFM (see above for details)
%       - 'onlyErrorFlag', [false] flag indicating if only the error is
%         needed when nFrame==2 && method==Inf && ~isProj
%       - 'nItrSba', [100] number of bundle adjustment iterations
%       - 'tolAffine', [1e-5] tolerance for affine upgrade as explained in
%                          affineUpgrade
%       - 'tolMetric', [1e-5] tolerance for affine upgrade as explained in
%                          metricUprgade
%
% OUTPUTS 1
%  anim      - Animation object
%
% OUTPUTS 2
%  err       - error, acoording to Gold Standard for Affine camera matrix
%
% EXAMPLE
%
% See also
%
% Vincent's Structure From Motion Toolbox      Version 3.1.1
% Copyright (C) 2008-2011 Vincent Rabaud.  [vrabaud-at-cs.ucsd.edu]
% Please email me if you find bugs, or have suggestions or questions!
% Licensed under the GPL [see external/gpl.txt]

[ isCalibrated doAffineUpgrade doMetricUpgrade K KFull method ...
  onlyErrorFlag nItrSBA tolAffine tolMetric ] = getPrmDflt( varargin, ...
  { 'isCalibrated' false 'doAffineUpgrade' false 'doMetricUpgrade' false ...
  'K' [] 'KFull' [] 'method' inf 'onlyErrorFlag' false 'nItrSBA' 100 ...
  'tolAffine' 1e-5 'tolMetric' 1e-5}, 1);

% Remove points that do not have two views with no NaN's
good_point_mask = sum(~any(isnan(W),1),3) >=2;
W = W(:,good_point_mask,:);

nFrame = size(W,3); nPoint = size(W,2);

if doMetricUpgrade; doAffineUpgrade=true; end
if nFrame==2
  doAffineUpgrade=false;
  if ~isCalibrated; doMetricUpgrade=false; end
end

if size(W,1)==3; W = normalizePoint(W,3); end

% create an animation object that wil contain the output
anim=Animation(); anim.isProj=isProj; WOri=W;

% If calibrated, apply inv(K)
if ~isempty(K); anim.K=K; end
if ~isempty(KFull); anim.KFull=KFull; end
KFull=anim.KFull;
if ~isempty(KFull)
  isCalibrated=true;
  anim.KFull=[];
  % unapply the internal parameter matrix to the measurements
  if size(KFull,3)==1
    W = normalizePoint(multiTimes(inv(KFull),normalizePoint(W,-3),1.2),3);
  else
    invKFull = zeros(3,3,nFrame);
    for i=1:nFrame; invKFull(:,:,i) = inv(KFull(:,:,i)); end
    W = normalizePoint(multiTimes(inv(KFull),normalizePoint(W,-3),2),3);
  end
end

if isProj
  [P,S] = computeSMFromWProjective( W, method, isCalibrated );
else
  [P,S] = computeSMFromWAffine( W, method, onlyErrorFlag );
end
if onlyErrorFlag; anim=P; return; end

% fill the Animation object with the results
anim.S=S; anim.P=P; anim.W=W;

% perform an affine upgrade if requested
H=[];
if anim.isProj
  % do bundle adjustment
  if nItrSBA > 0
    anim = bundleAdjustment( anim, 'nItr', nItrSBA );
  else
    anim = bundleAdjustment( anim, 'nItr', 20 );
  end
  
  % modify P so that P(:,:,1)==eye(3,4)
  anim=anim.setFirstPRtToId();
  
  % if we want to do an affine upgrade
  if doAffineUpgrade
    % only apply the quasi affine upgrade if we are under octave
    % (as Yalmip does not work there)
    if exist('OCTAVE_VERSION','builtin')==5
      warning('Only performing a quasi-affine upgrade under Octave');
      [ HEye, Hqa ] = affineUpgrade(anim);
      H=HEye*Hqa;
    else
      [ HEye, Hqa, pInf ] = affineUpgrade(anim, 'tol', tolAffine);
      % if the camera is calibrated or if we do not do a metric upgrade
      % stop here
      if ~doMetricUpgrade
        % apply H
        H=eye(4); H(4,1:3)=-pInf;
        H=HEye*H;
      end
    end
    if doMetricUpgrade
      % perform a metric upgrade if requested
      if exist('OCTAVE_VERSION','builtin')==5
        H=metricUpgrade(anim, 'isCalibrated', isCalibrated);
      else
        [ H, anim.KFull ]=metricUpgrade(anim, 'pInf', pInf, ...
          'isCalibrated', isCalibrated, 'tol', tolMetric);
      end
    end
  end
else
  if doMetricUpgrade
    H=metricUpgrade(anim, 'isCalibrated', isCalibrated, 'method', method);
  end
end

% apply the homography H
if ~isempty(H)
  anim.P=multiTimes(anim.P,H,1);
  anim.S=normalizePoint(H\normalizePoint(anim.S,-4),4);
end

% recover rotations and translations
if doMetricUpgrade && (exist('OCTAVE_VERSION','builtin')~=5 || isCalibrated)
  P=anim.P;
  if ~isCalibrated && ~isempty(anim.KFull)
    if size(anim.KFull,3)==1; P=multiTimes(inv(anim.KFull),P,1.2);
    else P=multiDiv(anim.KFull,P,2);
    end
  end
  R=zeros(3,3,nFrame);
  if anim.isProj
    for i=1:nFrame
      P(:,:,i)=P(:,:,i)/nthroot(det(P(:,1:3,i)),3);
      R(:,:,i)=rotationMatrix(P(:,1:3,i));
    end
    anim.t=reshape(P(:,4,:),3,nFrame);
  else
    for i=1:nFrame
      R(:,:,i)=rotationMatrix(rotationMatrix(P(1:2,1:3,i)));
    end
    anim.t=reshape(P(1:2,4,:),2,nFrame); anim.t(3,:)=0;
  end
  anim.R=R;
  
  anim=anim.setFirstPRtToId();
end

% update the mask if any
mask=reshape(any(~isnan(W),1),nPoint,nFrame);
if any(mask==0); anim.mask=mask; end

% do a final bundle adjustment
if nItrSBA > 0; anim = bundleAdjustment( anim, 'nItr', nItrSBA ); end

W = anim.W; S = anim.S;
anim.W = zeros(2, length(good_point_mask), nFrame);
anim.S = zeros(3, length(good_point_mask), nFrame);
anim.W = W; anim.S = S;

% Re-assign the KFull from the input arguments if any
if ~isempty(KFull); anim.KFull=KFull; anim.W=WOri; end
end
