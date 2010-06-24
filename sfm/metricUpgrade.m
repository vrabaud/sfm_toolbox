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
%       - 'pInf' [3 x 1] plane at inifinty found in affineUpgrade
%       - 'method' only used with an orthographic camera
%                   if 0, simple metric constraints
%                   if inf, metric constraints soolve mode exactly with SDP
%       - 'nItr' number of iterations in the branch and bound algorithm
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

[ isCalibrated method nItr ] = getPrmDflt( varargin, ...
  { 'isCalibrated' false 'method' inf 'nItr', 20}, 1);

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
  H=eye(4); H(4,1:3)=-pInf;
  return
end

% Yalmip does not work under octave sorry :(
if exist('OCTAVE_VERSION','builtin')==5; return; end

% Chandraker IJCV 2009
% globally optimal metric upgrade
HInfinity=anim.P(:,1:3,:)-multiTimes(anim.P(:,4,:), pInf',1);

% free variables for the convex/concave relaxations
omega=sdpvar(3,3);
omegaCol=[omega(1,1);omega(1,2);omega(1,3);omega(2,2);omega(2,3)];
omegaColOne=[omegaCol;1];
nu=sdpvar(3,3,nFrame);
nuCol=reshape([nu(1,1,:);nu(1,2,:);nu(1,3,:);nu(2,2,:);nu(2,3,:);...
  nu(3,3,:)],[6,nFrame]);
lam=sdpvar(1,nFrame); a=sdpvar(1,nFrame);

solverSetting = sdpsettings('solver','sedumi,sdpa,csdp,*','verbose',0, ...
  'cachesolvers',1);

FIni=[ omega >= 0, omega(3,3)==1 ];

% omega=KK' so we can deduce bounds on values of omega
kl=[0,-1000,0,-1000,-1000]; ku=[2000,1000,2000,1000,1000];
K=sdpvar(3,3);
F=[ kl(:)<=K([1:3,5:6,9])<=ku(:) ];
% l and u are the bounds of omega, where the values are
l=zeros(5,1); u=l;
for j=1:2
  for i=j:3
    solvesdp( F, K(j,:)*K(:,i), solverSetting );
    l(i+(j-1)*3)=double(K(j,:))*double(K(:,i));
    solvesdp( F, -K(j,:)*K(:,i), solverSetting );
    u(i+(j-1)*3)=-double(K(j,:))*double(K(:,i));
  end
end

% [ o1 o2 o3, o2 o4 o5; o3 o5 1 ]
omegaColBest=(l+u)/2; currBest=[Inf]; lowerBound=[Inf];
for itr=1:nItr
  % perform branch and bound
  [l,u,omegaColBest,currBest,lowerBound,newInd]=bnbBranch(l,u,omegaColBest,...
    currBest,lowerBound);

  % compute the lower bound/current best for each interval
  for ind=newInd
    F=FIni+[ l(:,ind)<=omegaCol<=u(:,ind) ];
    % compute Li and Ui
    L=zeros(1,nFrame); U=zeros(1,nFrame);
    for i=1:nFrame
      crit=HInfinity(3,:,i)*omega*HInfinity(3,:,i)';
      diagno = solvesdp( F, -crit, solverSetting );
      L(i)=1/double(crit);
      diagno = solvesdp( F, crit, solverSetting );
      U(i)=1/double(crit);
    end
    % add constraints
    F=FIni;
    aa=omegaColOne*U+[l(:,ind);1]*(lam-U);
    bb=omegaColOne*L+[u(:,ind);1]*(lam-L);
    cc=omegaColOne*L+[l(:,ind);1]*(lam-L);
    dd=omegaColOne*U+[u(:,ind);1]*(lam-U);
    F=F+[ cc(:)<=nuCol(:)<=aa(:), dd(:)<=nuCol(:)<=bb(:), L<=lam<=U ];

%      % define all the SOCP constraints
%      for i=1:nFrame
%        tmp=HInfinity(:,:,i)*nu(:,:,i)*HInfinity(:,:,i)';
%        F=F+cone(tmp(:)-omega(:),a(i));
%      end
%      crit=norm(a)^2;

    % define all the SOCP constraints
    crit=0;
    for i=1:nFrame
      tmpHInfinity(:,:,i)*nu(:,:,i)*HInfinity(:,:,i)';
      crit=crit+norm(tmp(:)-omega(:))^2;
    end

    % solve for omega
    diagno = solvesdp( F, crit, solverSetting );
    lowerBound(ind)=double(crit);

    % compute the current best
    omegaColBestTmp=double(omegaCol);
    currBestTmp=criterion(omegaColBestTmp,double(lam),HInfinity);
    if currBestTmp<currBest(ind)
      currBest(ind)=currBestTmp; omegaColBest(:,ind)=omegaColBestTmp;
    end
  end
  % remove intervals for which the lower bound is higher than the current
  % best of another intervals
  badInterval=find(lowerBound>min(currBest));
  lowerBound
  currBest
  [l;u]
  omegaColBest
  min(currBest)
  l(:,badInterval)=[]; u(:,badInterval)=[]; omegaColBest(:,badInterval)=[];
  lowerBound(:,badInterval)=[]; currBest(:,badInterval)=[];
  lowerBound
  currBest
  [l;u]
  omegaColBest
  min(currBest)
end
% figure out K from omega
[ disc, ind ]=min(currBest);
K=chol(triToFull([omegaColBest(:,ind);1]));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function res=g(a,b)
% perform some simple product for Kanade factorization
res=[ a(1)*b(1) a(1)*b(2)+a(2)*b(1) a(1)*b(3)+a(3)*b(1) a(2)*b(2) ...
  a(2)*b(3)+a(3)*b(2) a(3)*b(3) ];
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function res=triToFull(x)
% phaving a 6x1 or 1x6 vector, build the full 3x3 symmetric metrix
res=[x(1:3);x([2,4,5]);x([3,5,6])];
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function res=criterion(omegaCol,lam,HInfinity)
omega=triToFull([omegaCol;1]);
bestTmp=bsxfun(@minus, omega, bsxfun(@times,...
  reshape(lam,1,1,length(lam)), multiTimes(HInfinity,...
  multiTimes(omega,HInfinity,1.2),2)));
res=norm(bestTmp(:))^2;
end
