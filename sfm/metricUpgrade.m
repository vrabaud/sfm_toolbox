function [ H, K ] = metricUpgrade(anim, varargin)
% Perform an affine upgrade
%
% Given projection matrices and a 3D projective structure in anim,
% perform an affine upgrade
%
% USAGE
%   anim = affineUpgrade(anim, isCalibrated)
%
% INPUTS
%  anim          - Animation object with P and S filled
%  varargin   - list of paramaters in quotes alternating with their values
%       - 'isCalibrated' flag indicating if the camera is calibrated or not
%                  (K=eye(3))
%       - 'pInf' [3 x 1] plane at infinity found in affineUpgrade
%       - 'method' only used with an orthographic camera
%                   if 0, simple metric constraints
%                   if inf, metric constraints solve more exactly with SDP
%       - 'tol' [1e-5] difference between lower bound and optimal value
%               to stop at in the branch and bound algorithm
%
% OUTPUTS
%  H             - homography to apply to get projection matrices
%                  of the form K*[R|t] for projective cameras
%                  and [ K(1:2,:)*[R|t]; 0 0 0 1] for orthographic cameras
%  K             - calibration matrix
%
% EXAMPLE
%
% See also
%
% Vincent's Structure From Motion Toolbox      Version NEW
% Copyright (C) 2008-2010 Vincent Rabaud.  [vrabaud-at-cs.ucsd.edu]
% Please email me if you find bugs, or have suggestions or questions!
% Licensed under the GPL [see external/gpl.txt]

[ isCalibrated method tol pInf ] = getPrmDflt( varargin, ...
  { 'isCalibrated' false 'method' inf 'tol' 1e-5 'pInf' []}, 1);

P=anim.P; nFrame=anim.nFrame;
H=[]; K=[];

if ~anim.isProj
  if isCalibrated && (method==0 || exist('OCTAVE_VERSION','builtin')==5)
    % Tomasi Kanade with the metric constraint
    H=getPCloseToRotation(P,2);
  else
    if exist('OCTAVE_VERSION','builtin')==5
      error(['Cannot do an orthographic uncalibrated metric upgrade ' ...
        'under Octave']);
    end
    % solve for K and H such that P(1:2,1:3,i)*H=K(1:2,1:2)*R(1:2,:,i),
    % with R(:,:,i) a rotation and K is the calibration matrix of the form:
    % [ k1 k2 0; 0 k3 0; 0 0 1 ]
    % We therefore must have:
    % P(1:2,1:3,i)*H*H'*P(1:2,1:3,i)'=K(1:2,1:2)*K(1:2,1:2)'
    HHt=sdpvar(3,3);
    if isCalibrated && method==inf
      KKt=eye(2); K=eye(3); F=set(HHt>=0);
    else
      KKt=sdpvar(2,2); F=set(HHt>=0)+set(KK>=0);
    end
    % define all the SOCP constraints
    obj=sdpvar(2*nFrame,2);
    for i=1:nFrame
      obj(2*i-1:2*i,:)=P(1:2,1:3,i)*HHt*P(1:2,1:3,i)'-KKt;
    end
    diagno = solvesdp( F, norm(obj), sdpsettings('solver', ...
      'sedumi,sdpa,csdp,*','verbose',0) )';
    % compute the calibration matrix
    if ~isCalibrated; K=chol(double(KKt)); K(3,3)=1; end
    % compute H (up to a rotation ambiguity)
    H=chol(double(HHt)); H(4,4)=1;
  end
  return
end

% simple case where the camera is projective and calibrated
if isCalibrated
  if exist('OCTAVE_VERSION','builtin')~=5
    % we have pInf from an affine upgrade so it's all good
    H=eye(4); H(4,1:3)=-reshape(pInf,1,3);
    return;
  end
  % we do like Tomasi-Kanade with the metric constraint
  H=getPCloseToRotation(P,3);
  return;
end
if exist('OCTAVE_VERSION','builtin')==5
  error(['Cannot do an affine uncalibrated metric upgrade ' ...
    'under Octave']);
end

% solve for omega the conventional way, just to get bounds on it
% compute homographies at infinity
HInfinity=anim.P(:,1:3,:)-bsxfun(@times,anim.P(:,4,:),reshape(pInf,1,3,1));
for i=1:nFrame
  HInfinity(:,:,i)=HInfinity(:,:,i)/nthroot(det(HInfinity(:,:,i)),3);
end
% we are solving for sum norm(omega-Hinf*omega*Hinf')
A=zeros(9*nFrame,9);
for i=1:nFrame
  A(9*i-8:9*i,:)=eye(9)-kron(HInfinity(:,:,i),HInfinity(:,:,i));
end
% impose the symmetry on omega
A(:,2)=A(:,2)+A(:,4);
A(:,3)=A(:,3)+A(:,7);
A(:,6)=A(:,6)+A(:,8);
A(:,[4,7,8])=[];
[U,S,V]=svd(A); omega=[V(1:3,end)';V(2,end),V(4:5,end)';...
  V(3,end),V(5:6,end)'];

% impose semidefinitess artificially
[U,S,V]=svd(omega);

S(S<0)=0; omega=U*S*U';
omega=omega/omega(3,3);
% get the KK' transformation (kindof like cholesky which is R'R)
K=kFromOmega(omega);

% optimize K a bit more
KVect=[K(1,:)';K(2,2:3)'];
KVect=fminunc(@(x)criterionFromKVect(x,HInfinity),KVect,...
  optimset('GradObj','off','Hessian','off','Algorithm',...
  'active-set','Display','off'));

% deduce bounds on K (+-200%, totally empirical)
% K=[k1,k2,k3;0,k4,k5;0,0,1] and not K=[k1,k2,k4;0,k3,k5;0,0,1] as usual
kl=KVect;
ku=kl+2*abs(kl); kl=kl-2*abs(kl);
kl(1)=0; kl(4)=0;

% Chandraker IJCV 2009
% globally optimal metric upgrade
HInfinity=anim.P(:,1:3,:)-bsxfun(@times,anim.P(:,4,:),reshape(pInf,1,3,1));

% free variables for the convex/concave relaxations
omega=sdpvar(3,3);
omegaCol=[omega(1,1);omega(1,2);omega(1,3);omega(2,2);omega(2,3)];
omegaColOne=[omegaCol;1];
nu=sdpvar(3,3,nFrame);
nuCol=reshape([nu(1,1,:);nu(1,2,:);nu(1,3,:);nu(2,2,:);nu(2,3,:);...
  nu(3,3,:)],[6,nFrame]);
lam=sdpvar(1,nFrame); coneXMain=sdpvar(6*nFrame,1); coneYMain=sdpvar(1,1);

% build the objective function
for i=1:nFrame
  tmp=HInfinity(:,:,i)*nu(:,:,i)*HInfinity(:,:,i)';
  coneXMain(6*i-5:6*i)=[omega([1,5,9])'-tmp([1,5,9])';...
    sqrt(2)*(omega([4,7,8])'-tmp([4,7,8])')];
end

% build the constraints that we will respect for every interval
FIni=[ omega >= 0, omega(3,3)==1, cone(coneXMain,coneYMain) ];
solverSetting = sdpsettings('solver','sedumi,sdpa,csdp,*','verbose',0, ...
  'cachesolvers',1);

% [ o1 o2 o3, o2 o4 o5; o3 o5 1 ]
kVectBest=(kl+ku)/2; currBest=Inf; lowerBound=-Inf; isRefined=false;

[kl,ku,kVectBest,currBest,lowerBound,isRefined]=bnbRefine(...
  @(x)criterionFromKVect(x,HInfinity),kl,ku,kVectBest,currBest,...
  lowerBound,isRefined);

ticId = ticStatus('metric upgrade iterations',1,1);
nItr=100;
for itr=1:nItr
  % perform branch and bound
  [kl,ku,kVectBest,currBest,lowerBound,isRefined,newInd]=bnbBranch(...
    kl,ku,kVectBest,currBest,lowerBound,isRefined);
  
  % compute the lower bound/current best for each interval
  for ind=newInd
    % let us deduce L and U from the bounds on K
    % compute Li and Ui differently from the paper so that we don't have
    % to solve a bunch of SDP's. That's why we are doing bnb on k and not
    % omega
    L=zeros(1,nFrame); U=zeros(1,nFrame);
    kFullBound=zeros(3,3,2);
    kFullBound(:,:,1)=KVectToKFull(kl(:,ind));
    kFullBound(:,:,2)=KVectToKFull(ku(:,ind));
    for i=1:nFrame
      % we want to get the minimum/maximum of the elements of h3*K
      h3Bound=repmat(HInfinity(3,:,i),[1,1,2]);
      bound=reshape(productBound(h3Bound,kFullBound),3,2);
      % deduce the bounds on the norm
      bound=normBound(bound);
      % use 1e-6 so that the inverse is not too big
      L(i)=1/max([bound(2),1e-6]);
      U(i)=1/max([bound(1),1e-6]);
    end
    
    % omega=KK' so we can deduce bounds on values of omega
    [l,u]=omegaBound(kl(:,ind),ku(:,ind));
    F=FIni+[ l<=omegaCol<=u ];
    
    % add constraints
    aa=omegaColOne*U+[l;1]*(lam-U);
    bb=omegaColOne*L+[u;1]*(lam-L);
    cc=omegaColOne*L+[l;1]*(lam-L);
    dd=omegaColOne*U+[u;1]*(lam-U);
    F=F+[ cc(:)<=nuCol(:)<=aa(:), dd(:)<=nuCol(:)<=bb(:), L<=lam<=U ];
    
    % solve for omega
    diagno=solvesdp( F, coneYMain, solverSetting );
    lowerBound(ind)=double(coneYMain)^2;
    if ~isRefined(ind)
      KBestTmp=kFromOmega(double(omega));
      kVectBestTmp=[KBestTmp(1,:)';KBestTmp(2,2:3)'];
      % as bounds on omega are deduced from bounds on K, and as there is
      % no equivalence between the two, we need to make sure KVectBest
      % is the bounds for k
      kVectBest(:,ind)=min([max([kVectBestTmp,kl(:,ind)],[],2),...
        ku(:,ind)],[],2);
      currBest(ind)=criterionFromKVect(kVectBest(:,ind),HInfinity);
    end
  end
  % refine results if needed
  [kl,ku,kVectBest,currBest,lowerBound,isRefined]=bnbRefine(...
    @(x)criterionFromKVect(x,HInfinity),kl,ku,kVectBest,currBest,...
    lowerBound,isRefined);
  
  tocStatus( ticId, itr/nItr );
  if abs(min(currBest)-min(lowerBound))<tol; break; end
end
% figure out K from omega
[ disc, ind ]=min(currBest);
K=KVectToKFull(kVectBest(:,ind));
H=eye(4); H(1:3,1:3)=K; H(4,1:3)=-reshape(pInf,1,3)*K;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function H=getPCloseToRotation(P,nRow)
% Find H such that the top nRow rows of P*H is a bunch of rotations
% P is 3*4*nFrame and nRow is 2 or 3
% we therefore want P(1:nRow,:)*H*H'*P(1:nRow,:)'=eye(nRow), or again
% kron(P(1:nRow,:),P(1:nRow,:))*vect(H*H')=vect(eye(nRow))
nFrame=size(P,3); PKron=zeros(2*nRow,4*4,nFrame);
for row=1:nRow
  for i=1:4
    PKron(nRow*(row-1)+1:nRow*row,4*i-3:4*i,:)=bsxfun(@times,...
      P(row,i,:),P(1:nRow,:,:));
  end
end
A=eye(nRow); A=repmat(A(:),[nFrame,1]);
HHt=reshape(reshape(permute(PKron,[1,3,2]),nRow*nRow*nFrame,16)\A,4,4);
% now, artificially make sure it is symmetric
HHt=(HHt+HHt')/2;
% now, artificially make sure it is positive semi-definite
[ U S V ]=svd(HHt); S=diag(U'*HHt*U); SNeg=S<0; S(SNeg)=0;
H=U*diag(sqrt(S)); H(:,SNeg)=[]; H(:,4)=[0;0;0;1];
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function res=KVectToKFull(KVect)
% having the 5x1 KVect, build the full K matrix
res=[KVect(1:3)';0,KVect([4:5])';0,0,1];
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function K=kFromOmega(omega)
% make sure omega is positive semi-definite (should be very close, up to
% the SDP solver)
[U,S,V]=svd(omega); S(S<0)=0; omega2=U*S*V';
% perform homemade cholesky decomposition as the matrix can be semidefinite
L=zeros(3,3);
omega2=omega2(end:-1:1,end:-1:1);
for j=1:3
  for i=j:3
    if i==j
      L(i,i)=sqrt(max(omega2(i,i)-sum(L(i,1:i-1).^2),0));
    else
      L(i,j)=(omega2(i,j)-L(i,1:j-1)*L(j,1:j-1)')/L(j,j);
    end
  end
end
%  L=full(cholinc(sparse(omega(end:-1:1,end:-1:1)'),'0'));
% deduce K
L=L';
K=L(end:-1:1,end:-1:1)';
% just for numerical stability
K=K/K(3,3);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function omega=omegaFromKVect(KVect)
% compute omega from KVect [5x1]
K=zeros(3,3,size(KVect,2)); K(3,3,:)=1;
K(1,1,:)=KVect(1,:);
K(1,2,:)=KVect(2,:);
K(1,3,:)=KVect(3,:);
K(2,2,:)=KVect(4,:);
K(2,3,:)=KVect(5,:);
omega=multiTimes(K,K,2.2);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [l,u]=omegaBound(kl,ku)
% get bounds on omega, from bounds on K
bound=zeros(5,2);
klFull=KVectToKFull(kl); kuFull=KVectToKFull(ku);
% figure out the bounds on the diagonal elements
bound(1,:)=normBound([klFull(1,:)',kuFull(1,:)']);
bound(4,:)=normBound([klFull(2,:)',kuFull(2,:)']);

% figure out the bounds on the (1,2) element
tmp1=zeros(1,3,2); tmp2=zeros(3,1,2);
tmp1(1,:,1)=klFull(1,:); tmp1(1,:,2)=kuFull(1,:);
tmp2(:,1,1)=klFull(2,:)'; tmp2(:,1,2)=kuFull(2,:)';
bound(2,:)=reshape(productBound(tmp1,tmp2),1,2);

% the last two elements are easy to bound
bound([3,5],:)=[kl([3,5]), ku([3,5])];
l=bound(:,1); u=bound(:,2);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function bound=productBound(A,B)
% A is m x n x 2 and B is n x p x 2, bound is m x p x 2
tmp=sort(bsxfun(@times,reshape(A,size(A,1),size(A,2),1,2),...
  reshape(B,1,size(B,1),size(B,2),2)),4);
bound=reshape(sum(tmp,2),size(A,1),size(B,2),2);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function bound=normBound(A)
% A is m x 2, bound is 1 x 2
bound=sort(A.^2,2);
bound(A(:,1)<0 & A(:,2)>0,1)=0;
bound=sum(bound,1);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function res=criterionFromKVect(KVect,HInfinity)
% computes the criterion we are optimizing from KVect [5x1]
omega=omegaFromKVect(KVect);
% compute the optimal lambdas
nFrame=size(HInfinity,3); nSample=size(KVect,2);
% compute the criterion
res=zeros(1,nSample);
for i=1:nFrame
  HOmegaHt=multiTimes(HInfinity(:,:,i),multiTimes(omega,...
    HInfinity(:,:,i)',1),1.2);
  res=res+sum(reshape(omega-bsxfun(@rdivide,HOmegaHt,HOmegaHt(3,3,:)),...
    9,nSample).^2,1);
end
end
