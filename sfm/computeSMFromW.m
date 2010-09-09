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
% Returns the structure and camera parameters in an Animation object
%
% USAGE
%   anim = computeSMFromW( isProj, W, 'method',method )
%   anim = computeSMFromW( isProj, W,'method',Inf,'onlyErrorFlag',true )
%
% INPUTS
%  isProj     - flag indicating if the camera is projective
%  W          - [] [ 2 x nPoint x nFrame ] 2D projected features
%  varargin   - list of paramaters in quotes alternating with their values
%       - 'isCalibrated' [false] flag indicating if the cameras are
%                        calibrated (intrinsic parameters are identity)
%       - 'doAffineUpgrade' [false] flag indicating if the affine upgrade
%                        is computed
%       - 'doMetricUpgrade' [false] flag indicating if the metric upgrade
%                        is computed
%       - 'K',[] [3 x 1 ], [ 3 x nFrame ] calibration parameters
%                       (or [5 x 1 ], [ 5 x nFrame ] when projective)
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
% Vincent's Structure From Motion Toolbox      Version 3.0
% Copyright (C) 2008-2010 Vincent Rabaud.  [vrabaud-at-cs.ucsd.edu]
% Please email me if you find bugs, or have suggestions or questions!
% Licensed under the GPL [see external/gpl.txt]

[ isCalibrated doAffineUpgrade doMetricUpgrade K KFull method ...
  onlyErrorFlag nItrSBA tolAffine tolMetric ] = getPrmDflt( varargin, ...
  { 'isCalibrated' false 'doAffineUpgrade' false 'doMetricUpgrade' false ...
  'K' [] 'KFull' [] 'method' inf 'onlyErrorFlag' false 'nItrSBA' 100 ...
  'tolAffine' 1e-5 'tolMetric' 1e-5}, 1);

nFrame = size(W,3); nPoint = size(W,2);

if doMetricUpgrade; doAffineUpgrade=true; end
if nFrame==2
  doAffineUpgrade=false;
  if ~isCalibrated; doMetricUpgrade=false; end
end

P=zeros(3,4,nFrame); P(:,:,1)=eye(3,4); if ~isProj; P(3,3:4)=[0 1]; end

if size(W,1)==3; W = normalizePoint(W,3); end

% create an animation object that wil contain the output
anim=Animation(); anim.isProj=isProj; anim.W=W;

% If calibrated, apply inv(K)
if ~isempty(K); anim.K=K; end
if ~isempty(KFull); anim.KFull=KFull; end
if ~isempty(anim.K)
  isCalibrated=true;
  % unapply the calibration to the measurements
  KFull=anim.KFull;
  if size(KFull,3)==1
    invKFull = repmat( inv(KFull), [ 1 1 nFrame ] );
    KFull = repmat( KFull, [ 1 1 nFrame ] );
  else
    invKFull = zeros(3,3,nFrame);
    for i=1:nFrame; invKFull(:,:,i) = inv(KFull(:,:,i)); end
  end
  W = normalizePoint( multiTimes( invKFull, normalizePoint( W, -3 ), 2), 3 );
end

if isProj
  [P,S] = computeSMFromWProjective( W, method, isCalibrated );
else
  [P,S] = computeSMFromWProjective( W, method, onlyErrorFlag );
end

if onlyErrorFlag; anim=err; return; end

% fill the Animation object with the results
anim.S=S; anim.P=P;

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

  % if we want to do an affin upgrade, and if we are not in the case
  % of the essential matrix
  if doAffineUpgrade
    % only apply the quasi affine upgrade if we are under octave
    % (as Yalmip does not work there)
    if exist('OCTAVE_VERSION','builtin')==5
      [ HEye, Hqa ] = affineUpgrade(anim);
      H=HEye*Hqa;
    else
      [ HEye, Hqa, pInf ] = affineUpgrade(anim, 'tol', tolAffine);
      % if the camera is calibrated or if we do not do a metric upgrade
      % stop here
      if isCalibrated || ~doMetricUpgrade
        % apply H
        H=eye(4); H(4,1:3)=-pInf;
        H=HEye*H;
      else
        % perform a metric upgrade if requested
        if doMetricUpgrade
          [ H, anim.KFull ]=metricUpgrade(anim, 'pInf', pInf, ...
            'isCalibrated', isCalibrated, 'tol', tolMetric);
        end
      end
    end
  end
else
  if doMetricUpgrade
    if isCalibrated
      H=metricUpgrade(anim, 'isCalibrated', isCalibrated, 'method', method);
    elseif exist('OCTAVE_VERSION','builtin')~=5
      [ H, anim.KFull ]=metricUpgrade(anim, 'isCalibrated', ...
        isCalibrated, 'method', method);
    end
  end
end

% apply the homography H
if ~isempty(H)
  anim.P=multiTimes(anim.P,H,1);
  anim.S=normalizePoint(inv(H)*normalizePoint(anim.S,-4),4);
end

% recover rotations and translations
if doMetricUpgrade && (exist('OCTAVE_VERSION','builtin')~=5 || isCalibrated)
  P=anim.P; KFull=anim.KFull;
  if ~isCalibrated && ~isempty(KFull)
    if size(KFull,3)==1; P=multiTimes(inv(KFull),P,1.2);
    else P=multiDiv(KFull,P,2);
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
  
  % do a final bundle adjustment
  anim=anim.setFirstPRtToId();
  
  if nItrSBA > 0; anim = bundleAdjustment( anim, 'nItr', nItrSBA ); end
end
end
