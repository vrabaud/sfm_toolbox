function [ HEye, Hqa, pInf ] = affineUpgrade(anim, varargin)
% Perform an affine upgrade by recovering the plane at infinity
%
% Given projection matrices and a 3D projective structure in anim,
% perform an affine upgrade.
% Only quasi-affine is performed is performed on Octave though
%
% Basically an implementation of Chandraker IJCV 2009
%
% USAGE
%   anim = affineUpgrade(anim)
%
% INPUTS
%  anim          - Animation object with P and S filled
%  varargin   - list of paramaters in quotes alternating with their values
%     'doQuasiOnly'   - [false] flag indicating if we only want to do the quasi
%                       affine upgrade
%     'nItr'          - number of iterations in branch and bound
%
% OUTPUTS
%  HEye          - homography to apply to anim.P so that
%                  anim.P(:,:,1)*H=eye(3,4)
%  Hqa           - quasi affine homography after H has been applied
%  pInf          - best plane at infinity
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

[ doQuasiOnly nItr ] = getPrmDflt( varargin, ...
  { 'doQuasiOnly' false 'nItr' 20}, 1);

if exist('OCTAVE_VERSION','builtin')==5; doQuasiOnly=true; end

% make sure the first projection matrix is eye(3,4)
[anim,HEye]=anim.setFirstPRtToId();

PClean=anim.P; SHomogClean=normalizePoint(anim.S,-4); W=anim.W;
nPoint=anim.nPoint; nFrame=anim.nFrame;
Hqa=eye(4);

% notations from Chandraker IJCV 2009
% compute Hq that takes the projective frame to some quasi-affine frame
% (HZ2, p527)
WHat=W; WHat(3,:,:)=1;
w=squeeze(mean(multiTimes(PClean,SHomogClean,1)./WHat,1));
% change the signs of P's and X's
hasChanged=true;
while hasChanged
  hasChanged=false;
  if sum(sign(w(:,1)),1)<0
    SHomogClean=-SHomogClean; w=-w; hasChanged=true;
  end
  ind=sum(sign(w),1)<0;
  if nnz(ind) > 0
    PClean(:,:,ind)=-PClean(:,:,ind); w(:,ind)=-w(:,ind); hasChanged=true;
  end
  
  ind=sum(sign(w),2)<0;
  if nnz(ind) > 0
    SHomogClean(:,ind)=-SHomogClean(:,ind); w(ind,:)=-w(ind,:); hasChanged=true;
  end
end

% build the chirality inequalities
C=zeros(4,nFrame);
for j=1:nFrame
  for i=1:4; C(i,j)=(-1)^i*det(PClean(:,[1:i-1,i+1:4],j)); end
end

% remove points for which one w is negative
SHomogClean(:,any(w<0,2)) = [];

% find a solution to the chiral inequalities
for delta=-1:2:1
  % we want to minimize -d, with the parameters [v1,v2,v3,v4,d]
  A=-[SHomogClean';delta*C']; A(:,5)=1;
  if exist('OCTAVE_VERSION','builtin')==5
    [sol,fval] = linprog([0,0,0,0,-1]',A,zeros(size(A,1),1),zeros(0,5),...
      zeros(0,5),[-1,-1,-1,-1,-Inf]', [1,1,1,1,Inf]');
    if ~any(isna(sol)) && sol(5)>0; sol=solTmp; break; end
  else
    [solTmp,fval,exitflag] = linprog([0,0,0,0,-1]',A,zeros(size(A,1),1),...
      [],[],[-1,-1,-1,-1,-Inf]', [1,1,1,1,Inf]', [], ...
      optimset('Display','off'));
    if exitflag==1 && solTmp(5)>0; sol=solTmp; break; end
  end
end

% define Hq
Hq=eye(4); Hq(4,:)=sol(1:4)'; Hq(1,1)=delta/Hq(4,4);

% define Ha (HZ2, p528, algorithm 21.2)
% mean(Hq*X,2)=0 and svd(Hq*X)=1
HqS=normalizePoint(Hq*SHomogClean,4);
c=mean(HqS,2);
HqSCentered=bsxfun(@minus,HqS,c);
% make sure the moments are 1
[U,S]=eig(HqSCentered*HqSCentered');
S=diag(1./sqrt(diag(S)));
K=S*U'; if det(K) < 0; K = -K; end
Ha=[ K -K*c; 0 0 0 1 ];
% define the quasi-affine affinity that also centers the data
Hqa=Ha*Hq;

% Yalmip does not work under octave sorry :(
if doQuasiOnly; return; end

% Chandraker IJCV 2009

% figure out the bounds on vi (equation 8)
PQuasi=multiTimes(PClean,inv(Hqa),1);
SQuasi=normalizePoint(Hqa*SHomogClean,4);

% build the chirality inequalities
% algorithm 21.2 p528
C=zeros(4,nFrame);
for j=1:nFrame
  PQuasi(:,:,j)=PQuasi(:,:,j)/nthroot(det(PQuasi(:,1:3,j)),3);
  for i=1:4; C(i,j)=(-1)^i*det(PQuasi(:,[1:i-1,i+1:4],j)); end
end

A=-[C(1:3,:)';SQuasi']; b=[C(4,:)';ones(size(SQuasi,2),1)];
l=zeros(3,1); u=l;
for i=1:3
  for j=-1:2:1
    f=zeros(3,1); f(i)=j;
    if exist('OCTAVE_VERSION','builtin')==5
      [sol,fval] = linprog(f,A,b);
      if ~any(isna(sol)); exitflag=1; end
    else
      [sol,fval,exitflag] = linprog(f,A,b,...
        [],[],[],[],[],optimset('Display','off'));
    end
    if exitflag==1; if j<0; u(i)=sol(i); else l(i)=sol(i); end; end
  end
end
vBest=(u+l)/2;

% compute the greek coefficient
% the code below explains the chirality constraints
%syms lam a1 a2 a3 a4 a5 a6 a7 a8 a9 a10 a11 a12 p1 p2 p3;
%A=[a1 a2 a3 a4; a5 a6 a7 a8; a9 a10 a11 a12];
%alpha = -[ a4, a8, a12, - a11 - a6 - a1 ]*[p1;p2;p3;1];
%beta = [a2.*a8 + a12.*a3 - a11.*a4  - a4.*a6, a12.*a7 - a1.*a8 + ...
%  a4.*a5 - a11.*a8, - a1.*a12 - a12.*a6 + a10.*a8 + a4.*a9, a1.*a11 + ...
%  a1.*a6 - a2.*a5 + a11.*a6 - a10.*a7 - a3.*a9 ]*[p1;p2;p3;1];
%gamma = -[ a12.*a2.*a7 - a12.*a3.*a6 - a11.*a2.*a8 - a10.*a4.*a7 + ...
%  a11.*a4.*a6 + a10.*a3.*a8, -a1.*a12.*a7 + a12.*a3.*a5 + a1.*a11.*a8 - ...
%  a11.*a4.*a5 - a3.*a8.*a9 + a4.*a7.*a9, a1.*a12.*a6 - a12.*a2.*a5 - ...
%  a1.*a10.*a8 + a10.*a4.*a5 + a2.*a8.*a9 - a4.*a6.*a9, a11.*a2.*a5 - ...
%  a1.*a11.*a6 + a1.*a10.*a7 - a10.*a3.*a5 - a2.*a7.*a9 + a3.*a6.*a9 ]*...
%  [p1;p2;p3;1];
%simplify(lam^3 - alpha*lam^2 + beta*lam - gamma-det(lam*eye(3)-...
%  (A(:,1:3)-A(:,4)*[p1,p2,p3])))
P=PClean;
a1=reshape(P(1,1,:),nFrame,1); a2=reshape(P(1,2,:),nFrame,1);
a3=reshape(P(1,3,:),nFrame,1); a4=reshape(P(1,4,:),nFrame,1);
a5=reshape(P(2,1,:),nFrame,1); a6=reshape(P(2,2,:),nFrame,1);
a7=reshape(P(2,3,:),nFrame,1); a8=reshape(P(2,4,:),nFrame,1);
a9=reshape(P(3,1,:),nFrame,1); a10=reshape(P(3,2,:),nFrame,1);
a11=reshape(P(3,3,:),nFrame,1); a12=reshape(P(3,4,:),nFrame,1);
alpha = -[ a4, a8, a12, - a11 - a6 - a1 ];
beta = [ a2.*a8 + a12.*a3 - a11.*a4  - a4.*a6, a12.*a7 - a1.*a8 + ...
  a4.*a5 - a11.*a8, - a1.*a12 - a12.*a6 + a10.*a8 + a4.*a9, a1.*a11 + ...
  a1.*a6 - a2.*a5 + a11.*a6 - a10.*a7 - a3.*a9 ];
gamma = -[ a12.*a2.*a7 - a12.*a3.*a6 - a11.*a2.*a8 - a10.*a4.*a7 + ...
  a11.*a4.*a6 + a10.*a3.*a8, -a1.*a12.*a7 + a12.*a3.*a5 + a1.*a11.*a8 - ...
  a11.*a4.*a5 - a3.*a8.*a9 + a4.*a7.*a9, a1.*a12.*a6 - a12.*a2.*a5 - ...
  a1.*a10.*a8 + a10.*a4.*a5 + a2.*a8.*a9 - a4.*a6.*a9, a11.*a2.*a5 - ...
  a1.*a11.*a6 + a1.*a10.*a7 - a10.*a3.*a5 - a2.*a7.*a9 + a3.*a6.*a9 ];

% compute the product of the greek letters with Hqa
% greek will designate alpha,bet,gamma or [0,0,0,1]
greekHqat=zeros(nFrame,4,3);
for j=1:4
  switch j,
    case 1,
      greek=alpha;
    case 2,
      greek=beta;
    case 3,
      greek=gamma;
    case 4,
      greek=zeros(nFrame,4); greek(:,4)=1;
  end
  greekHqat(:,:,j)=greek*Hqa';
end

% free variables for the convex/concave relaxations
r=sdpvar(1,1); e=sdpvar(1,1); fg=sdpvar(nFrame,2);
v=sdpvar(3,1); % we'll set the 4th value to 1 by hand
t=sdpvar(nFrame,2); t1=sdpvar(nFrame,2,2);

% rotated Lorentz cone constraint
FIni=set(cone([2*(fg(:,1)-fg(:,2));r-e],r+e)) + ...
  set(e(:)>=0) + set(r(:)>=0);

currBest=Inf; lowerBound=Inf;
eps=1e-4;
solverSetting = sdpsettings('solver','sedumi,sdpa,csdp,*','verbose',0, ...
  'cachesolvers',0,'debug',1);
for itr=1:nItr
  % perform branch and bound
  [l,u,vBest,currBest,lowerBound,newInd]=bnbBranch(l,u,vBest,...
    currBest,lowerBound);
  
  % compute the lower bound/current best for each new interval
  for ind=newInd
    % add boundary conditions
    F=FIni+set(l(:,ind)<=v<=u(:,ind));
    
    % compute the bounds on a,b,c,d
    abcdBound=zeros(nFrame,2,4);
    for j=1:4
      % figure out the bounds of aj, bj, cj, d
      % get the smallest and biggest elements of greekHqat*v
      allProduct=sort(bsxfun(@times,greekHqat(:,:,j),reshape(...
        [l(:,ind) u(:,ind); 1 1],1,4,2)),3);
      abcdBound(:,:,j)=reshape(sum(allProduct,2),nFrame,2);
    end
    
    % add relaxation constraints on f and g
    for k=1:2
      F = F + relaxx13y(fg(:,k),greekHqat(:,:,k+2)*[v;1],...
        abcdBound(:,1,k+2),abcdBound(:,2,k+2),greekHqat(:,:,k)*[v;1],...
        abcdBound(:,1,k),abcdBound(:,2,k),t(:,k),t1(:,:,k));
    end
    % add the constraint on e
    F = F + concx83(e,greekHqat(1,:,4)*[v;1],abcdBound(1,1,4),...
      abcdBound(1,2,4));
    
    % solve the problem
    % yop, sad but that needs to be done
    clear mexsdpa;
    diagno = solvesdp(F,r,solverSetting);
    lowerBound(ind)=double(r);
    
    % make sure the best value is in the interval
    vBestTmp=min([max([double(v),l(:,ind)],[],2),u(:,ind)],[],2);
    
    % check for a few random values in the interval and check if they are
    % better than the value at the optimal v
    allV=[ vBestTmp, bsxfun(@plus,bsxfun(@times,rand(3,10000),...
      u(:,ind)-l(:,ind)), l(:,ind)) ];
    allMin=criterion(allV,greekHqat);
    [ disc, bestInd ]=min(allMin);
    vBestTmp=allV(:,bestInd);
    
    % perform gradient descent from the best v to find a better solution
    vBestTmp=fmincon(@(x)criterion(x,greekHqat),...
      vBestTmp,[],[],[],[],l(:,ind),u(:,ind),[],optimset(...
      'GradObj','off','Hessian','off','Algorithm','active-set'));
    % make sure the best value is in the interval
    vBestTmp=min([max([vBestTmp,l(:,ind)],[],2),u(:,ind)],[],2);
    % keep that v if the criterion is better than the current one
    currBestTmp=criterion(vBestTmp,greekHqat);
    if currBestTmp<currBest(ind)
      currBest(ind)=currBestTmp; vBest(:,ind)=vBestTmp;
    end
  end
  % remove intervals for which the lower bound is higher than the current
  % best of another interval
  mini=min(currBest);
  badInterval=find(lowerBound>min(currBest));

% the following is necessary because of roundups
  if ~isempty(badInterval) && length(badInterval)~=length(currBest)
    l(:,badInterval)=[]; u(:,badInterval)=[]; vBest(:,badInterval)=[];
    currBest(badInterval)=[]; lowerBound(badInterval)=[];
  end
  lowerBound
  currBest
  [l;u]
  vBest
  min(currBest)
  if size(l,2)==1 && norm(l-u,'fro')<eps; break; end
end
% figure out the plane at infinity
[disc,ind]=min(currBest);
pInf=Hqa'*[vBest(:,ind);1];
pInf=pInf(1:3)/pInf(4);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function conc=concx83(z,x,xl,xu)
% compute the concave relaxation for x^(8/3)
conc=[ z <= nthroot(xl,3).^8+(x-xl)./(xu-xl).*(nthroot(xu,3).^8-...
  nthroot(xl,3).^8) ];
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function constr=relaxxy(z,x,xl,xu,y,yl,yu)
% compute convex/concave relaxations for x*y
% precompute some data (because yalmip is slow ...)
xly=xl.*y; ylx=yl.*x; xuy=xu.*y; yux=yu.*x;
% compute the concave relaxation for x*y
constr=[ z <= xuy + ylx - xu.*yl; z <= xly + yux - xl.*yu ];
% compute the convex relaxation for x*y
constr=[ constr; z >= xly + ylx - xl.*yl; z >= xuy + yux - xu.*yu ];
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [a,b]=lineTangentCoeff(x)
% basically the line in equation A.8 for a point x
a=1/3*nthroot(x,3).^(-2); b=2/3*nthroot(x,3);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function constr=convx13yxposCase1(z,x,xl,xu,y,yl,yu)
% compute the convex relaxation for x^(1/3)*y when xl > 0 and yl > 0
% precompute lambda
a=1./(xu-xl); lam=a.*x-a.*xl;
% precompute the linear coefficients
a=nthroot(xl,3)-nthroot(xu,3); b=nthroot(xu,3).*y;
% add the constraints
constr=[ z >= a.*yu - lam.*(a.*yu) + b, z >= a.*y - lam.*(a.*yl) + b ];
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function constr=convx13yxposCase2(z,x,xl,xu,y,yl,yu)
% compute the convex relaxation for x^(1/3)*y when xl > 0 and yl < 0 < yu
% precompute some values as Yalmip can be slow
a1=(nthroot(xl,3)-nthroot(xu,3)).*yu./(xl-xu);
b1=(nthroot(xu,3).*xl-nthroot(xl,3).*xu).*yu./(xl-xu);
lam=(y-yl)./(yu-yl);
[a3,b3]=lineTangentCoeff(xl); [a4,b4]=lineTangentCoeff(xu);
% add constraints
constr=[];
for t1=[ (1-lam).*xu, x-lam.*xl ]
  for t2=[ a3.*t1 + (1-lam).*b3, a4.*t1 + (1-lam).*b4 ]
    constr=constr+[ z>=t2.*yl + a1.*(x-t1)+lam.*b1 ];
  end
end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function constr=convx13yxposCase3(z,x,xl,xu,y,yl,yu,t1)
% compute the convex relaxation for x^(1/3)*y when xl > 0 and yu < 0
lam=(y-yl)./(yu-yl);
% precompute lam.*xl and lam.*xu (because yalmip is slow)
lamxl=lam.*xl; lamxu=lam.*xu;
% constraints on t1
constr=[ xl-lamxl<=t1<=xu-lamxu; lamxl<=x-t1<=lamxu ];
% precompute values for constraints involving t2 and t3
[a3,b3]=lineTangentCoeff(xl); [a4,b4]=lineTangentCoeff(xu);
ylt2Set=[(yl.*a3).*t1+(1-lam).*(yl.*b3),(yl.*a4).*t1+(1-lam).*(yl.*b4)];
yut3Set=[(yu.*a3).*(x-t1)+lam.*(yu.*b3),(yu.*a4).*(x-t1)+lam.*(yu.*b4)];
% create the constraints
for ylt2=ylt2Set
  for yut3=yut3Set
    constr=constr+[ z>=ylt2+yut3 ];
  end
end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function constr=convx13yxpos(z,x,xl,xu,y,yl,yu,t1)
% compute the convex relaxation for x^(1/3)*y and xl>0
% case 1
ind=(0<yl);
if nnz(ind) >= 1
  constr=convx13yxposCase1(z(ind),x(ind),xl(ind),xu(ind),...
    y(ind),yl(ind),yu(ind));
else
  constr=[];
end
% case 2
ind=(yl<0) & (0<yu);
if nnz(ind) >= 1
  constr=constr+convx13yxposCase2(z(ind),x(ind),xl(ind),xu(ind),...
    y(ind),yl(ind),yu(ind));
end
% case 3
ind=(yu<0);
if nnz(ind) >= 1
  constr=constr+convx13yxposCase3(z(ind),x(ind),xl(ind),xu(ind),...
    y(ind),yl(ind),yu(ind),t1(ind));
end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function constr=convx13yCaseI(z,x,xl,xu,y,yl,yu,t1)
% convex relaxation for x^(1/3)*y, when xl>0 or xu<0 (Case I of the paper)
% case where xl<0
ind=(xl>0);
if nnz(ind) >= 1
  constr=convx13yxpos(z(ind),x(ind),xl(ind),xu(ind),y(ind),...
    yl(ind),yu(ind),t1(ind));
else
  constr=[];
end
% case where xu<0. Simply negate x and y
ind=(xu<0);
if nnz(ind) >= 1
  constr=constr+convx13yxpos(z(ind),-x(ind),-xu(ind),-xl(ind),-y(ind),...
    -yu(ind),-yl(ind),t1(ind));
end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function constr=convx13CaseII(x,xl,xu,t)
% convex relaxation for x^(1/3)*y, when xl<0<xu (Case II of the paper)
% only for t=x^(1/3) though
% t convex constraints
ind=(xl<=0) & (0<=xu) & (-xu/8>xl);
if nnz(ind)>=1
  [a,b]=lineTangentCoeff(xl(ind));
  %    constr=[ t(ind)>=(nthroot(xu(ind),-3).^2-2/3*...
  %      nthroot(xl(ind),3)./xu(ind)).*x(ind)+2/3*nthroot(xl(ind),3); ...
  %      t(ind)>=a.*x(ind)+b ];
  % the following violates Cauchy conditions but is bounded by the
  % commented one above that does, so we are safe for the convergence
  constr=[ t(ind)>=(nthroot(xu(ind),3).^(-2)*4/3.*(x(ind)-xu(ind))+ ...
    nthroot(xu(ind),3)); t(ind)>=a.*x(ind)+b ];
else
  constr=[];
end
% t convex constraints bis
ind=(xl<=0) & (0<=xu) & (-xu/8<=xl);
if nnz(ind)>=1
  % equation A.10
  a=(nthroot(xu(ind),3)-nthroot(xl(ind),3))./(xu(ind)-xl(ind));
  b=(xu(ind).*nthroot(xl(ind),3)-xl(ind).*nthroot(xu(ind),3))./...
    (xu(ind)-xl(ind));
  constr=[ constr; t(ind) >= a.*x(ind)+b ];
end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function constr=relaxx13y(z,x,xl,xu,y,yl,yu,t,t1)
% Case I: xl>0 or xu<0
% add the convex and the concave relaxation for x^(1/3)*y
constr=convx13yCaseI(z,x,xl,xu,y,yl,yu,t1(:,1)) + ...
  convx13yCaseI(-z,x,xl,xu,-y,-yu,-yl,t1(:,2));

% Case II: xl<0<xu
% convex and concave constraints on t
constr=constr + convx13CaseII(x,xl,xu,t) + convx13CaseII(-x,-xu,-xl,-t);

% add the constraints on z=ty
ind=(xl<=0) & (0<=xu);
if nnz(ind)>=1
  constr=[ constr; relaxxy(z(ind),t(ind),nthroot(xl(ind),3),...
    nthroot(xu(ind),3),y(ind),yl(ind),yu(ind)) ];
end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [res,grad]=criterion(v,greekHqat)
v(4,:)=1;
a=greekHqat(:,:,1)*v; b=greekHqat(:,:,2)*v;
c=greekHqat(:,:,3)*v; d=greekHqat(:,:,4)*v;
res=sum((nthroot(c,3).*a-nthroot(d,3).*b).^2,1)./nthroot(d(1,:),3).^8;
grad=[];
if nargout==2;
  % compute the gradient too
  deriv1=(bsxfun(@times,1/3*nthroot(c,3).^(-2).*a,greekHqat(:,1:3,3))+...
    bsxfun(@times,nthroot(c,3),greekHqat(:,1:3,1))-...
    bsxfun(@times,1/3*nthroot(d,3).^(-2).*b,greekHqat(:,1:3,4))-...
    bsxfun(@times,nthroot(d,3),greekHqat(:,1:3,2)));
  deriv1=sum(bsxfun(@times,deriv1,2*(nthroot(c,3).*a-nthroot(d,3).*b)./...
    nthroot(d(1,:),3).^8),1);
  grad=deriv1+(-8/3)*res*(greekHqat(1,1:3,4)/d(1));
  grad=grad';
end
end
