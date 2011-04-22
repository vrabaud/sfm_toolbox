function [P,S] = computeSMFromWProjective( W, method, isCalibrated )
% Compute SfM from measurements for projective camera
%
% This function should not be called directly: it does not contain
% pre-processing steps (removal of calibration matrix if known) or
% post-procssing steps (affine/metric upgrades) and bundle adjustment
% computeSMFromW should always be used instead.
%
% Essential matrix decomposition
% Reference: HZ2, p259, Result 9.19, and p. 294
% nFrame==2 && isCalibrated
% fast
%
% Normalized 8-point algorithm
% Reference: HZ2, p279 and alg 11.1 p. 282
% nFrame==2 && ~isCalibrated
% fast
%
% Projective Factorization
% Reference: HZ2, p445, Algorithm 18.2
% Sturm and Triggs ECCV 96
% A Factorization Based Algorithm for Multi-Image Projective Structure and
% Motion
% nFrame>2 && method==0
%
% Projective Factorization
% Reference: Iterative Extensions of the Sturm/Triggs Algorithm:
% Convergence and Nonconvergence, from Oliensis, Hartley, PAMI 07
% nFrame>2 && method==Inf
% slower
%
% If there are any missing entries:
% Low-Rank Matrix Fitting Based on Subspace Perturbation Analysis
% with Applications to Structure from Motion
% Hongjun Jia, Aleix M. Martinez, PAMI 08
%
% Returns the structure and camera parameters in an Animation object
%
% USAGE
%   anim = computeSMFromWProjective( W, method )
%
% INPUTS
%  W            - [ 2 x nPoint x nFrame ] 2D projected features
%  method       - [inf] method for performing SFM (see above for details)
%  isCalibrated - flag indicating if the camera is calibrated (only used
%                 for 2 frames)
%
% OUTPUTS
%  P          - [ 3 x 4 x nFrame ] projection matrices
%  S          - [ 3 x nPoint ] structure matrix
%
% EXAMPLE
%
% See also
%
% Vincent's Structure From Motion Toolbox      Version NEW
% Copyright (C) 2008-2010 Vincent Rabaud.  [vrabaud-at-cs.ucsd.edu]
% Please email me if you find bugs, or have suggestions or questions!
% Licensed under the GPL [see external/gpl.txt]

nFrame = size(W,3); nPoint = size(W,2);

P=zeros(3,4,nFrame); P(:,:,1)=eye(3,4);

WIsnan=isnan(W);
hasAnyNan=any(WIsnan(:));

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
if nFrame>2 && ismember(method,[ 0 Inf ]) && ~hasAnyNan
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

% Low-Rank Matrix Fitting Based on Subspace Perturbation Analysis
% with Applications to Structure from Motion
% Hongjun Jia, Aleix M. Martinez, PAMI 08
if nFrame>2 && ismember(method,[ 0 Inf ]) && hasAnyNan
  % Normalize coordinates
  [ q TW ] = normalizePoint(W, -Inf);
  HHat=reshape(permute(q,[1,3,2]),3*nFrame,nPoint);
  TInv=1./permute(sum(q.*q,1),[3,2,1]);
  Sk=HHat;
  for k = 1 : 30
    % Sk is HHat (homogeneous W) but scaled by lambda, initialized to 1
    [PStack,S]=lowRankDecomposition(Sk,4);
    
    % make sure the columns of PStack are orthonormal
    [PStack,SDiag,V]=svd(PStack,'econ');
    S = SDiag*V'*S;
    
    P = permute(reshape(PStack,3,nFrame,4),[1,3,2]);
    %sumTot=0;
    for j=1:nPoint
      mask=~isnan(TInv(:,j));
      Cj=permute(sum(bsxfun(@times,q(:,j,mask),P(:,:,mask)),1),[3,2,1]);
      % Solve the maximum of the Rayleigh quotient
      % it is the highest eigenvector of A*v=lam*(1./TjInv)*v
      % or heighest eigen value of TjInv*A*v=lam*v
      A=Cj*Cj';
      [lamjTmp,val]=eigs(bsxfun(@times,TInv(mask,j),A),1,'lm');
      lamj=zeros(1,1,nFrame); lamj(mask)=lamjTmp;
      %sumTot=sumTot+lamjTmp'*Cj*Cj'*lamjTmp/(lamjTmp'*diag(1./ ...
      % TInv(mask,j))*lamjTmp);
      Sk(:,j)=reshape(bsxfun(@times,q(:,j,:),lamj),3*nFrame,1);
    end
    %sumTot;
  end
  [PStack,S]=lowRankDecomposition(Sk,4);
  
  PTmp = PStack;
  S = normalizePoint(S,4);
  P = zeros(3,4,nFrame);
  for i=1:nFrame; P(:,:,i) = TW(:,:,i)\PTmp(3*i-2:3*i,:); end
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
