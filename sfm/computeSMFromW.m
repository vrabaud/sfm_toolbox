function anim = computeSMFromW( isProj, W, varargin )
% Recover rigid structure and motion from 2D calibrated measurements W
%
% Normalized 8-point algorithm
% Reference: HZ2, p279 and alg 11.1 p. 282
% nFrame==2 && method==0 && isProj
%
% Gold Standard for Affine camera matrix
% Reference: HZ2, p351, Algorithm 14.1
% nFrame==2 && ~isProj && onlyErrorFlag
%
% Tomasi Kanade without the metric constraint
% Affine camera matrix, MLE estimation (Tomasi Kanade)
% Reference: HZ2, p437, Algorithm 18.1
% nFrame>=2 && method==0 && ~isProj
%
% Tomasi Kanade with the metric constraint
% nFrame>=2 && method==0 && ~isProj && isCalibrated
%
% Projective Factorization
% Reference: HZ2, p445, Algorithm 18.2
% Sturm and Triggs ECCV 96
% A Factorization Based Algorithm for Multi-Image Projective Structure and
% Motion
% nPoint>=2 && method==0 && isProj && ~isCalibrated
%
% Projective Factorization
% Reference: Iterative Extensions of the Sturm/Triggs Algorithm:
% Convergence and Nonconvergence, from Oliensis, Hartley, PAMI 07
% nPoint>=2 && method==Inf && isProj && ~isCalibrated
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
%       - 'doAffineUpgrade' [true] flag indicating if the affine upgrade
%                        is computed
%       - 'doMetricUpgrade' [true] flag indicating if the metric upgrade
%                        is computed
%       - 'K',[] [3 x 3 ] or [ 3 x 3 x nFrame ] calibration matrices
%       - 'method', [inf] method for performing SFM (see above for details)
%       - 'onlyErrorFlag', [false] flag indicating if only the error is
%         needed when nFrame==2 && method==Inf && ~isProj
%       - 'nItrSBA', [100] number of bundle adjustment iterations
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

[ isCalibrated doAffineUpgrade doMetricUpgrade K method onlyErrorFlag ...
  nItrSBA ] = getPrmDflt( varargin, { 'isCalibrated' false ...
  'doAffineUpgrade' true 'doMetricUpgrade' true 'K' [] 'method' inf ...
  'onlyErrorFlag' false 'nItrSBA' 100 }, 1 );

nFrame = size(W,3); nPoint = size(W,2);

P=zeros(3,4,nFrame); P(:,:,1)=eye(3,4); if ~isProj; P(3,3:4)=[0 1]; end

if size(W,1)==3
  W = normalizePoint(W,3);
  W = W( 1:2, :, : );
end
WOri = W;

% If calibrated, apply inv(K)
if ~isempty(K)
  isCalibrated = true;
  if size(K,3)==1
    invK = repmat( inv(K), [ 1 1 nFrame ] );
    K = repmat( K, [ 1 1 nFrame ] );
  else
    invK = zeros(3,3,nFrame);
    for i=1:nFrame; invK(:,:,i) = inv(K(:,:,i)); end
  end
  W = normalizePoint( multiTimes( invK, normalizePoint( W, -3 ), 2), 3 );
end

if isProj
  % Normalized 8-point algorithm
  % Reference: HZ2, p279 and alg 11.1 p. 282
  if nFrame==2 && method==0
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
    else
      % Essential matrix decomposition
      % Reference: HZ2, p259, Result 9.19, and p. 294
      [ U disc V ] = svd(F);
      
      % Check which of the 4 possibilities gives a point in front of both
      % cameras
      WW = [ 0 -1 0; 1 0 0; 0 0 1 ];
      maxCount = 0;
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
        if det(P(:,1:3,2))<0; P(:,1:3,2) = -P(:,1:3,2); end
        
        C = zeros(3,2); v = C;
        for j = 1 : 2
          % camera center(Reference: HZ2, p158-161
          [ disc disc VV ] = svd(P(:,:,j)); C(:,j) =VV(1:3,end)'/VV(4,end);
          % principal axis
          M = P(:,1:3,j); v(:,j) = det(M)*M(3,:)';
        end
        % compute the image of a point
        S = computeSFromWM( true, W, P, 'method', 0);
        maxCountTmp = sum((v(:,1)'*(S-repmat(C(:,1),1,nPoint)))>0) + ...
          sum((v(:,2)'*(S-repmat(C(:,2),1,nPoint)))>0);
        if maxCountTmp>maxCount; P2 = P(:,:,2); maxCount = maxCountTmp; end
      end
      P(:,:,2) = P2;
    end
    
    S = computeSFromWM( true, W, P, 'method', 0 );
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
      [ Wkm1, Wkm1Hat ] = projSturmTriggsWkWkHat(lam,W,nFrame,nPoint);
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
      lam=sum(WWkm1Hat,1)./WSquaredSum;
      
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

% create the Animation object
anim=Animation(); anim.isProj=isProj; anim.W=WOri; anim.S=S; anim.P=P;
if ~isempty(K); anim.K = K; end

% do bundle adjustment
if nItrSBA > 0; anim = bundleAdjustment( anim, 'nItr', nItrSBA ); end

% perform an affine upgrade if requested
if doAffineUpgrade; anim = affineUpgrade(anim, isCalibrated); end

% perform a metric upgrade if requested
if doMetricUpgrade; anim = metricUpgrade(anim, isCalibrated); end

% do a final bundle adjustment
%if nItrSBA > 0; anim = bundleAdjustment( anim, 'nItr', nItrSBA ); end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [ Wk, WkHat ] =projSturmTriggsWkWkHat(lam,W,nFrame,nPoint)
% perform the computation of Wk and its rank 4 approximation
% in Storm Triggs
Wk = bsxfun(@times,lam,W);
[ U S V ] = svd( reshape(permute(Wk,[1,3,2]),3*nFrame,nPoint),...
  'econ' );
WkHat = permute(reshape(U(:,1:4)*S(1:4,1:4)*V(:,1:4)',3,nFrame,...
  nPoint),[1,3,2]);
end
