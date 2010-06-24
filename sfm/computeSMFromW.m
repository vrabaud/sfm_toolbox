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
%       - 'nItrAff', [20] number of iterations in the affine upgrade
%       - 'nItrMetr', [20] number of iterations in the metric upgrade
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
  onlyErrorFlag nItrSBA nItrAff nItrMetr ] = getPrmDflt( varargin, ...
  { 'isCalibrated' false 'doAffineUpgrade' false 'doMetricUpgrade' false ...
  'K' [] 'KFull' [] 'method' inf 'onlyErrorFlag' false 'nItrSBA' 100 ...
  'nItrAff' 20 'nItrMetr' 20}, 1);

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
  % Normalized 8-point algorithm
  % Reference: HZ2, p279 and alg 11.1 p. 282
  if nFrame==2
    % Normalize input data
    [x T]=normalizePoint(W(:,:,1),Inf);
    [xp Tp]=normalizePoint(W(:,:,2),Inf);
    
    A=[xp([1 1],:).*x; xp(1,:); xp([2 2],:).*x;xp(2,:);x;ones(1,nPoint)]';
    
    [U,S,V]=svd(A,0); F=reshape(V(:,end),[3,3])';
    [U,S,V]=svd(F,0);
    F=U*diag([S(1,1) S(2,2) 0])*V';
    
    F=[ Tp; 0 0 1 ]'*F*[ T; 0 0 1 ];
    
    if ~isCalibrated
      P(:,:,2) = convertPF([],F,true);
      S = computeSFromWM( true, W, P, 'method', 0 );
    else
      % Essential matrix decomposition
      % Reference: HZ2, p259, Result 9.19, and p. 294
      [ U disc V ] = svd(F);
      if det(U)<0; U=-U; end; if det(V)<0; V=-V; end
      
      % Check which of the 4 possibilities gives a point in front of both
      % cameras
      WW = [ 0 -1 0; 1 0 0; 0 0 1 ];
      maxCount=0;
      for i = 1 : 4
        switch i
          case 1,
            P(:,:,2) = [ U*WW*V' U(:,3) ];
          case 2,
            P(:,:,2) = [ U*WW*V' -U(:,3) ];
          case 3,
            P(:,:,2) = [ U*WW'*V' U(:,3) ];
          case 4,
            P(:,:,2) = [ U*WW'*V' -U(:,3) ];
        end
        
        C = zeros(3,2); v = C;
        for j = 1 : 2
          % camera center (Reference: HZ2, p158-161)
          M = P(:,1:3,j);
          C(:,j)=-M\P(:,4,j);
          % principal axis
          v(:,j) = det(M)*M(3,:)';
        end
        % compute the image of a point
        STmp = computeSFromWM( true, W, P, 'method', 0);
        maxCountTmp = sum((v(:,1)'*bsxfun(@minus,STmp,C(:,1)))>0) + ...
          sum((v(:,2)'*bsxfun(@minus,STmp,C(:,2)))>0);
        if maxCountTmp>maxCount
          P2 = P(:,:,2); maxCount = maxCountTmp; S=STmp;
        end
      end
      P(:,:,2) = P2;
    end
  end
  
  % Projective Factorization
  % Reference: Iterative Extensions of the Sturm/Triggs Algorithm:
  % Convergence and Nonconvergence, from Oliensis, Hartley, PAMI 07
  if nFrame>2 && ismember(method,[ 0 Inf ])
    % Normalize coordinates
    [ W T ] = normalizePoint(W,-Inf);
    
    % Initialize
    lam = ones( 1, nPoint, nFrame );
    WSquaredSum = sum(W.^2,1);
    C0Const = sum(WSquaredSum(:));
    for n=1:50
      % Stage 1
      if method==0
        for k = 1 : 10
          % normalize each column and then each row
          lam=bsxfun(@rdivide, lam, sqrt(sum(lam.^2,2)));
          lam=bsxfun(@rdivide, lam, sqrt(sum(lam.^2,3)));
        end
      end
      
      % get the best rank 4 approximation
      % Wkm1 is for Wk minus 1
      [Wkm1,Wkm1Hat]=projSturmTriggsWkm1Wkm1Hat(lam,W);
      WWkm1Hat = sum(W.*Wkm1Hat,1);

      if method==Inf
        if n==1
          mu = norm(Wkm1(:)-Wkm1Hat(:))^2/norm(Wkm1(:))^4*1.1;
        else
          C0 = mu*C0Const;
          C1 = mu*sum(WWkm1Hat(:));
          C2 = WWkm1Hat.^2./WSquaredSum; C2 = mu*sum(C2(:));
          C3 = mu*sum(Wkm1Hat(:).^2);
        end
      end
      
      % Stage 2
      % get the optimal lambdas
      lam=WWkm1Hat./WSquaredSum;
      
      % Stage 3
      if method==Inf && n>1
        a = roots( [ C0, -(C0^2-2*C1), -(2*C0*C3-C2), ...
          -(4*C1*C3-2*C2*C0), C0*C3^2-2*C2*C3, ...
          2*C1*C3^2-C2^2, C2*C3^2 ] );
        a = real(a(abs(imag(a))<1e-10)); a = a(a>=0);
        
        if length(a)>1
          zkm1=C3*C0/C2;
          if abs(zkm1-1)>eps
            a = a( sign(zkm1-1)==sign(a/sqrt(C3)-1) );
          else a = sqrt(C2/C0)+1;
          end
        end
        lam = ( a + lam )*sqrt(a/(a^2*C0+2*a*C1+C2));
        %        [a,sqrt(a/(a^2*C0+2*a*C1+C2)), mean(lam(:)), std(lam(:))]
      end
    end
    
    % De-normalize coordinates
    Wk = bsxfun(@times,lam,W);
    [ U S V ] = svd( reshape(permute(Wk,[1,3,2]),3*nFrame,nPoint),...
      'econ' );
    P = U(:,1:4)*S(1:4,1:4); S = V(:,1:4)';
    
    PTmp = P;
    S = normalizePoint(S,4);
    P = zeros(3,4,nFrame);
    for i=1:nFrame; P(:,:,i) = T(:,:,i)\PTmp(3*i-2:3*i,:); end
  end
else
  % Gold Standard for Affine camera matrix
  % Reference: HZ2, p351, Algorithm 14.1
  % Just the error is computed
  if nFrame==2 && onlyErrorFlag
    % Normalize input data
    A=[ W(:,:,2); W(:,:,1) ]';
    A=bsxfun(@minus,A,mean(A,1));
    
    N = solveLeastSqAx( A, [], 1 );
    
    errFrame = N'*A';
    err = norm(errFrame)^2;
  else
    % Tomasi Kanade without/with the metric constraint
    % Affine camera matrix, MLE estimation (Tomasi Kanade)
    % Reference: HZ2, p437, Algorithm 18.1
    if nFrame>=2
      WStack = reshape( permute( W, [ 1 3 2 ] ), [], nPoint );
      
      t = mean( WStack, 2 );
      WStack = bsxfun(@minus, WStack, t );
      [ U S V ] = svd( WStack );
      P = permute(reshape(U(:,1:3),2,nFrame,3),[1 3 2]);
      P(3,4,:) = 1; P(1:2,4,:) = reshape( t, 2, 1, nFrame );
      
      S=bsxfun(@times,[ S(1,1); S(2,2); S(3,3) ], V(:,1:3)');
    end
  end
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
      [ HEye, Hqa, pInf ] = affineUpgrade(anim, 'nItr', nItrAff);
      % if the camera is calibrated or if we do not do a metric upgrade
      % stop here
      if isCalibrated || ~doMetricUpgrade
        % apply H
        H=eye(4); H(4,1:3)=-pInf;
        H=HEye*H;
      else
        % perform a metric upgrade if requested
        if doMetricUpgrade
          [ anim.KFull, H ]=metricUpgrade(anim, pInf, 'isCalibrated', ...
            isCalibrated, 'nItr', nItrMetr);
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [ Wkm1, Wkm1Hat ]=projSturmTriggsWkm1Wkm1Hat(lam,W)
% perform the computation of Wk and its rank 4 approximation
% in Storm Triggs
nFrame=size(W,3); nPoint=size(W,2);
Wkm1 = bsxfun(@times,lam,W);
[ U S V ] = svd( reshape(permute(Wkm1,[1,3,2]),3*nFrame,nPoint), 'econ' );
Wkm1Hat = permute(reshape(U(:,1:4)*S(1:4,1:4)*V(:,1:4)',3,nFrame,...
  nPoint),[1,3,2]);
end
