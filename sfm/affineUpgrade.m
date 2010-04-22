function anim = affineUpgrade(anim)
% Perform an affine upgrade
%
% Given projection matrices and a 3D projective structure in anim,
% perform an affine upgrade
%
% USAGE
%   anim = affineUpgrade(anim, isCalibrated)
%
% INPUTS
%  anim       - Animation object with P and S filled
%
% OUTPUTS
%  anim      - Animation object such that 
%
% EXAMPLE
%
% See also
%
% Vincent's Structure From Motion Toolbox      Version 3.0
% Copyright (C) 2008-2010 Vincent Rabaud.  [vrabaud-at-cs.ucsd.edu]
% Please email me if you find bugs, or have suggestions or questions!
% Licensed under the GPL [see external/gpl.txt]

if ~anim.isProj; return; end

P=anim.P; S=anim.S; W=anim.W; nPoint=anim.nPoint; nFrame=anim.nFrame;
SOne=S; SOne(4,:)=1;

% Chandraker IJCV 2009
% compute Hq that takes the projective frame to some quasi-affine frame
% (HZ2, p527)
WHat=W; WHat(3,:,:)=1;
w=squeeze(mean(multiTimes(P,[S;ones(1,nPoint)],1)./WHat,1));
% change the signs of P's and X's
for j=2:nFrame
  if sum(sign(w(:,j))-sign(w(:,1)))>0.5*nPoint
    P(:,:,j)=-P(:,:,j); w(:,j)=-w(:,j);
  end
end
for i=1:nPoint
  if nnz(w(i,:)<0) > nnz(w(i,:)>0); S(:,i)=-S(:,i); w(i,:)=-w(i,:); end
end
% build the chirality inequalities
C=zeros(4,nFrame);
for j=1:nFrame
  for i=1:4; C(i,j)=(-1)^i*det(P(:,[1:i-1,i+1:4],j)); end
end
% find a solution to the chiral inequalities
for delta=[-1:2:1]
  % we want t minimize -d, with the parameters [v1,v2,v3,v4,d]
  A=[SOne';C']; A(:,5)=0;
  [sol,fval,exitflag] = linprog([0,0,0,0,-1],-A,0,[],[],...
    [-1,-1,-1,-1,-Inf], [1,1,1,1,Inf]);
  if exitflag==1; break; end
end
% define Hq
Hq=eye(4); Hq(4,:)=sol; Hq(1,1)=delta/Hq(4,4);

% define Ha
Sq=Hq*SOne; Sq=bsxfun(@div,Sq(1,:3,:),Sq(4,:));
Ha=eye(4); Ha(1:3,4)=-mean(Sq,2);
Hqa=Ha*Hq;

% 
%  syms lam a1 a2 a3 a4 a5 a6 a7 a8 a9 a10 a11 a12 p1 p2 p3;
%  A=[a1 a2 a3 a4; a5 a6 a7 a8; a9 a10 a11 a12];
%  alpha = -(a4*p1 + a8*p2 + a12*p3 - a11 - a6 - a1);
%  beta = (a2*a8 + a12*a3 - a11*a4  - a4*a6)*p1 + (a12*a7 - a1*a8 + a4*a5 - a11*a8)*p2 + (- a1*a12 - a12*a6 + a10*a8 + a4*a9)*p3 + a1*a11 + a1*a6 - a2*a5 + a11*a6 - a10*a7 - a3*a9;
%  gamma = - ((a12*a2*a7 - a12*a3*a6 - a11*a2*a8 - a10*a4*a7 + a11*a4*a6 + a10*a3*a8)*p1 + (-a1*a12*a7 + a12*a3*a5 + a1*a11*a8 - a11*a4*a5 - a3*a8*a9 + a4*a7*a9)*p2 + (a1*a12*a6 - a12*a2*a5 - a1*a10*a8 + a10*a4*a5 + a2*a8*a9 - a4*a6*a9)*p3 + a11*a2*a5 - a1*a11*a6 + a1*a10*a7 - a10*a3*a5 - a2*a7*a9 + a3*a6*a9);
%  simplify(lam^3 - alpha*lam^2 + beta*lam - gamma-det(lam*eye(3)-(A(:,1:3)-A(:,4)*[p1,p2,p3])))
