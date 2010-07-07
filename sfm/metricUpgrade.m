function [ H, K ] = metricUpgrade(anim, pInf, varargin)
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
%                   if inf, metric constraints soolve mode exactly with SDP
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
% Vincent's Structure From Motion Toolbox      Version 3.0
% Copyright (C) 2008-2010 Vincent Rabaud.  [vrabaud-at-cs.ucsd.edu]
% Please email me if you find bugs, or have suggestions or questions!
% Licensed under the GPL [see external/gpl.txt]

[ isCalibrated method tol ] = getPrmDflt( varargin, ...
  { 'isCalibrated' false 'method' inf 'tol', 1e-5}, 1);

P=anim.P; S=anim.S; W=anim.W; nPoint=anim.nPoint; nFrame=anim.nFrame;
H=[]; K=[];

if ~anim.isProj
  if isCalibrated && (method==0 || exist('OCTAVE_VERSION','builtin'))
    % Tomasi Kanade with the metric constraint
    G=zeros( 3*nFrame, 6 );
    for i=1:nFrame
      G(i,:) = g( P(1,:,i), P(1,:,i) );
      G(nFrame+i,:) = g( P(2,:,i), P(2,:,i) );
      G(2*nFrame+i,:) = g( P(1,:,i), P(2,:,i) );
    end
    % Solve for the square root matrix
    l = G\[ ones(2*nFrame,1); zeros(nFrame,1) ];
    L = [ l(1:3)'; l(2) l(4:5)'; l(3) l(5:6)' ];
    [ U S V ] = svd(L); S(S<0) = 0; H = U*sqrt(S)*V';

    H(4,4)=1;
  else
    % solve for K and A such that P(1:2,1:3,i)*H=K(1:2,1:2)*R(1:2,:,i),
    % with R(:,:,i) a rotation and K is the calibration matrix of the form:
    % [ k1 k2 0; 0 k3 0; 0 0 1 ]
    % We therefore must have:
    % P(1:2,1:3,i)*H*H'*P(1:2,1:3,i)'=K(1:2,1:2)*K(1:2,1:2)'
    HH=sdpvar(3,3); a=sdpvar(1,nFrame);
    if isCalibrated && method==inf
      KK=eye(2); K=eye(3); F=set(HH>=0);
    else
      KK=sdpvar(2,2); F=set(HH>=0)+set(KK>=0);
    end
    % define all the SOCP constraints
    for i=1:nFrame
      F=F+cone(P(1:2,1:3,i)*HH*P(1:2,1:3,i)'-KK,a(i));
    end
    diagno = solvesdp( F, sum(a), sdpsettings('solver', ...
        'sdpa,csdp,sedumi,*','verbose',0) );
    % compute the calibration matrix
    if ~isCalibrated; K=chol(double(KK)); K(3,3)=1; end
    % compute H (up to a rotation ambiguity)
    H=chol(double(HH)); H(4,4)=1;
  end
  return
end

% simple case where the camera is projective and calibrated
if anim.isProj && isCalibrated
  H=eye(4); H(4,1:3)=-reshape(pInf,1,3);
  return
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
  A(9*i-8:9*i,:)=kron(eye(3),eye(3))-...
    kron(HInfinity(:,:,i),HInfinity(:,:,i));
end
[U,S,V]=svd(A); omega=reshape(V(:,end),3,3); omegaOri=omega;
% impose semidefinitess artificially
[U,S,V]=svd(omega); S(S<0)=0; omega=U*S*U';
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
if kl(1)<0; kl(1)=0; end; if kl(4)<0; kl(4)=0; end

% Yalmip does not work under octave sorry :(
if exist('OCTAVE_VERSION','builtin')==5; return; end

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
lam=sdpvar(1,nFrame);

FIni=[ omega >= 0, omega(3,3)==1 ];
solverSetting = sdpsettings('solver','sedumi,sdpa,csdp,*','verbose',0, ...
  'cachesolvers',1);

% [ o1 o2 o3, o2 o4 o5; o3 o5 1 ]
kVectBest=(kl+ku)/2; currBest=[Inf]; lowerBound=[Inf];

ticId = ticStatus('metric upgrade iterations',1,1);
nItr=100;
for itr=1:nItr
  % perform branch and bound
  [kl,ku,kVectBest,currBest,lowerBound,newInd]=bnbBranch(kl,ku,...
    kVectBest,currBest,lowerBound);

  % compute the lower bound/current best for each interval
  for ind=newInd
    % let us deduce L and U from the bounds on K
    % compute Li and Ui differently from the paper so that we don't have
    % to solve a bunch of SDP's. That's why we are doing bnb on k and not
    % omega
    L=zeros(1,nFrame); U=zeros(1,nFrame);
    klFullT=KVectToKFull(kl(:,ind))'; kuFullT=KVectToKFull(ku(:,ind))';
    for i=1:nFrame
      h3=HInfinity(3,:,i)';
      KtH3=zeros(3,2);
      for j=1:3
        % we want to get the minimum/maximum of the elements of K'*h
        tmp=[klFullT(j,:).*h3';kuFullT(j,:).*h3'];
        tmp=sort(tmp,1);
        KtH3(j,:)=sum(tmp,2)';
        % now deduce the bounds on the squared elements of K'*h
        if KtH3(j,1)<=0 && KtH3(j,2)>=0
          KtH3(j,:)=[0,max(KtH3(j,:).^2)];
        else
          KtH3(j,:)=sort(KtH3(j,:).^2);
        end
      end
      % deduce the bounds on the norm
      % use 1e-6 for roundup errors
      L(i)=1/double(max([sum(KtH3(:,2)),1e-6]));
      U(i)=1/double(max([sum(KtH3(:,1)),1e-6]));
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

    % define all the SOCP constraints
    crit=0;
    for i=1:nFrame
      tmp=HInfinity(:,:,i)*nu(:,:,i)*HInfinity(:,:,i)';
      crit=crit+norm(tmp(:)-omega(:))^2;
    end

    % solve for omega
    diagno=solvesdp( F, crit, solverSetting );
    lowerBound(ind)=double(crit);
    KBestTmp=kFromOmega(double(omega));
    KVectBestTmp=[KBestTmp(1,:)';KBestTmp(2,2:3)'];
    
    % refine the best value
    [KVectBestTmp,currBestTmp]=bnbRefine(kl(:,ind),ku(:,ind),...
      KVectBestTmp,@(x)criterionFromKVect(x,HInfinity));

    if currBestTmp<currBest(ind)
      currBest(ind)=currBestTmp; kVectBest(:,ind)=KVectBestTmp;
    end
  end
  % remove intervals for which the lower bound is higher than the current
  % best of another intervals
  badInterval=find(lowerBound>min(currBest));
  if ~isempty(badInterval) && length(badInterval)~=length(currBest)
    kl(:,badInterval)=[]; ku(:,badInterval)=[]; kVectBest(:,badInterval)=[];
    lowerBound(:,badInterval)=[]; currBest(:,badInterval)=[];
  end
%    lowerBound
%    currBest
%    min(currBest)
  tocStatus( ticId, itr/nItr );
  if abs(min(currBest)-min(lowerBound))<tol; break; end
end
% figure out K from omega
[ disc, ind ]=min(currBest);
K=KVectToKFull(kVectBest(:,ind));
H=eye(4); H(1:3,1:3)=K; H(4,1:3)=-reshape(pInf,1,3)*K;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function res=g(a,b)
% perform some simple product for Kanade factorization
res=[ a(1)*b(1) a(1)*b(2)+a(2)*b(1) a(1)*b(3)+a(3)*b(1) a(2)*b(2) ...
  a(2)*b(3)+a(3)*b(2) a(3)*b(3) ];
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function res=KVectToKFull(KVect)
% having the 5x1 KVect, build the full K matrix
res=[KVect(1:3)';0,KVect([4:5])';0,0,1];
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function K=kFromOmega(omega)
% extract K from omega
K=full(cholinc(sparse(omega(end:-1:1,end:-1:1)'),'inf'));
K=K(end:-1:1,end:-1:1)';
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
l=zeros(5,1); u=l;
ind=1;
klFull=KVectToKFull(kl); kuFull=KVectToKFull(ku);
% figure out the bounds on the diagonal elements
for i=[1,2]
  tmp=[ klFull(i,:).*klFull(i,:); kuFull(i,:).*kuFull(i,:) ];
  tmp(3,:)=tmp(1,:).*~(klFull(i,:)<0 & 0<kuFull(i,:));
  if i==1; j=1; else; j=4; end
  l(j)=sum(min(tmp,[],1)); u(j)=sum(max(tmp,[],1));
end
tmp=[ klFull(1,:).*klFull(2,:); kuFull(1,:).*klFull(2,:); ...
  kuFull(1,:).*kuFull(2,:) ];
l(i)=sum(min(tmp,[],1)); u(i)=sum(max(tmp,[],1));
% the last two elements are easy to bound
l([3,5])=kl([3,5]); u([3,5])=ku([3,5]);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function res=criterion(omegaCol,lam,HInfinity)
% computes the criterion we are optimizing using omega and lambda
omega=[omegaCol(1:3);omegaCol([2,4,5]);omegaCol([3,5]),1];
bestTmp=bsxfun(@minus, omega, bsxfun(@times,...
  reshape(lam,1,1,length(lam)), multiTimes(HInfinity,...
  multiTimes(omega,permute(HInfinity,[2,1,3]),1.2),2)));
res=norm(bestTmp(:))^2;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function res=criterionFromKVect(KVect,HInfinity)
% computes the criterion we are optimizing from KVect [5x1]
omega=omegaFromKVect(KVect);
% compute the optimal lambdas
nFrame=size(HInfinity,3); nSample=size(KVect,2);
lam=zeros(nFrame,nSample);
for i=1:nFrame
  lam(i,:)=reshape(1./multiTimes(multiTimes(HInfinity(3,:,i),omega,1.2),...
    HInfinity(3,:,i)',1),1,nSample);
end
% compute the criterion
res=zeros(1,nSample);
for i=1:nFrame
  resTmp=omega - bsxfun(@times,...
    reshape(lam(i,:),1,1,nSample), multiTimes(HInfinity(:,:,i),...
    multiTimes(omega,HInfinity(:,:,i)',1),1.2));
  res=res+reshape(sum(sum(resTmp.^2,1),2),1,nSample);
end
end
