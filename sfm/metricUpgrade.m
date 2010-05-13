function [ H, K ] = metricUpgrade(anim, isCalibrated, vBest)
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
%  isCalibrated  - flag indicating if the camera is calibrated or not
%                  (K=eye(3))
%  vBest         - optimal v found in affineUpgrade
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

P=anim.P; S=anim.S; W=anim.W; nPoint=anim.nPoint; nFrame=anim.nFrame;
K=[];

if ~anim.isProj
  if isCalibrated
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
    HH=sdpvar(3,3); KK=sdpvar(2,2); a=sdpvar(1,nFrame);
    F=set(HH>=0)+set(KK>=0);
    % define all the SOCP constraints
    for i=1:nFrame
      F=F+set(cone(P(1:2,1:3,i)*HH*P(1:2,1:3,i)-KK,a(i)));
    end
    diagno = solvesdp( F, sum(a), sdpsettings('solver', ...
        'sdpa,csdp,sedumi,*','debug',2) );
    % compute the calibration matrix
    K=chol(double(KK)); K(3,3)=1;
    % compute H (up to a rotation ambiguity)
    H=chol(double(HH)); H(4,4)=1;
  end
  return
end

% simple case where the camera is projective and calibrated
if anim.isProj && isCalibrated
  piInf=Hqa'*vBest; H=eye(4); H(4,1:3)=-piInf;
  return
end
return
% Chandraker IJCV 2009
% globally optimal metric upgrade
HInfinity=anim.P(:,1:3,:)-multiTimes(anim.P(:,4,:),[v(1),v(2),v(3),1],1);

% free variables for the convex/concave relaxations
omega=sdpvar(3,3); omegaCol=[omega(1,1);omega(1,2);omega(1,3);...
  omega(2,2);omega(2,3)];
nu=sdpvar(5,nFrame); lam=sdpvar(1,nFrame);
FIni=set(omega >= 0) + set(omega(3,3)==1);

l=repmat(-100,5,1); u=-l;
% l and u are the bounds of omega, where the values are
% [ o1 o2 o3, o2 o4 o5; o3 o5 1 ]
for nItr=1:10
  % perform branch and bound
  lNew=zeros(3,0); uNew=zeros(3,0);
  for i=1:size(l,2)
    % find the dimension to split onto
    [ disc, dim ]=max(u-l);
    uTmp=u(:,i); uTmp(dim)=(u(dim,i)+l(dim,i))/2;
    lNew(:,end+1)=l(:,i); uNew(:,end+1)=uTmp;
    lNew(:,end+1)=uTmp; uNew(:,end+1)=u(:,i);
  end
  l=lNew; u=uNew;

  % compute the lower bound/current best for each interval
  omegaBest=zeros(3,3,size(l,2));
  currBest=zeros(1,size(l,2));
  lowerBound=zeros(1,size(l,2));
  for ind=1:size(l,2)
    % compute Li and Ui
    L=zeros(1,nFrame); U=zeros(1,nFrame);
    F=FIni+set(l(:,ind)<=omegaCol<=u(:,ind));
    for i=1:nFrame
      crit=HInfinity(3,:,i)*omega*HInfinity(3,:,i)';
      diagno = solvesdp( F, crit, sdpsettings('solver', ...
        'sdpa,csdp,sedumi,*','debug',2) );
      L(i)=double(crit);
      diagno = solvesdp( F, -crit, sdpsettings('solver', ...
        'sdpa,csdp,sedumi,*','debug',2) );
      U(i)=-double(crit);
    end
    % add constraints
    F=F+set(nu<=omegaCol*U+l(:,ind)*lam-l(:,ind)*U);
    F=F+set(nu<=omegaCol*L+u(:,ind)*lam-u(:,ind)*L);
    F=F+set(nu>=omegaCol*L+l(:,ind)*lam-l(:,ind)*L);
    F=F+set(nu>=omegaCol*U+u(:,ind)*lam-u(:,ind)*U);
    F=F+set(l(:,ind)<=omegaCol<=u(:,ind));
    F=F+set(L<=lam<=U);

    % solve for omega
    crit=bsxfun(@minus, omega, multiTimes(HInfinity,...
      multiTimes(vi,HInfinity,2),2));
    diagno = solvesdp( F, crit, sdpsettings('solver', ...
        'sdpa,csdp,sedumi,*','debug',2) );
    omegaBest(:,:,ind)=[omegaCol(1:3); omegaCol([2,4,5]); ...
      omegaCol([3,5]), 1];
    
    % compute the current best
    bestTmp=bsxfun(@minus, double(omega), bsxfun(@times,...
      reshape(double(lam),1,1,nFrame), multiTimes(HInfinity,...
      multiTimes(double(omega),HInfinity,1.2),2)));
    curBest(i)=norm(bestTmp(:))^2;
    lowerBestTmp=double(crit);
    lowerBound(i)=norm(lowerBestTmp(:))^2;
  end
  % remove intervals for which the lower bound is higher than the current
  % best of another intervals
  badInterval=find(lowerBound>min(curBest));
  l(:,badInterval)=[]; u(:,badInterval)=[]; omegaBest(:,:,badInterval)=[];
end
% figure out K from omega
K=chol(omegaBest(:,:,1));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function res=g(a,b)
% perform some simple product for Kanade factorization
res=[ a(1)*b(1) a(1)*b(2)+a(2)*b(1) a(1)*b(3)+a(3)*b(1) a(2)*b(2) ...
  a(2)*b(3)+a(3)*b(2) a(3)*b(3) ];
end
