function [ Hqa, vBest ] = affineUpgrade(anim)
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
%
% OUTPUTS
%  Hqa           - quasi affine homography
%  vBest         - best v (to compute plane at infinity)
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
Hqa=eye(3); vBest=zeros(3,1);

% notations from Chandraker IJCV 2009
% compute Hq that takes the projective frame to some quasi-affine frame
% (HZ2, p527)
WHat=W; WHat(3,:,:)=1; SOne=[S;ones(1,size(S,2))];
w=squeeze(mean(multiTimes(P,SOne,1)./WHat,1));
% change the signs of P's and X's
hasChanged=true;
while hasChanged
  hasChanged=false;
  ind=sum(sign(w),1)<0;
  if nnz(ind) > 0
    P(:,:,ind)=-P(:,:,ind); w(:,ind)=-w(:,ind); hasChanged=true;
  end

  ind=sum(sign(w),2)<0;
  if nnz(ind) > 0
    SOne(:,ind)=-SOne(:,ind); w(ind,:)=-w(ind,:); hasChanged=true;
  end
end

% build the chirality inequalities
C=zeros(4,nFrame);
for j=1:nFrame
  for i=1:4; C(i,j)=(-1)^i*det(P(:,[1:i-1,i+1:4],j)); end
end

% remove points for which one w is negative
SOneClean=SOne; SOneClean(:,any(w<0,2)) = [];

% find a solution to the chiral inequalities
for delta=-1:2:1
  % we want to minimize -d, with the parameters [v1,v2,v3,v4,d]
  A=-[SOneClean';delta*C']; A(:,5)=1;
  if exist('OCTAVE_VERSION','builtin')==5
    [sol,fval] = linprog([0,0,0,0,-1]',A,zeros(size(A,1),1),zeros(0,5),...
      zeros(0,5),[-1,-1,-1,-1,-Inf]', [1,1,1,1,Inf]');
    if ~any(isna(sol)) && sol(5)>0; break; end
  else
    [sol,fval,exitflag] = linprog([0,0,0,0,-1]',A,zeros(size(A,1),1),...
      [],[],[-1,-1,-1,-1,-Inf]', [1,1,1,1,Inf]', [], ...
      optimset('Display','off'));
    if exitflag==1 && sol(5)>0; break; end
  end
end
% define Hq
Hq=eye(4); Hq(4,:)=sol(1:4)'; Hq(1,1)=delta/Hq(4,4);

% define Ha (HZ2, p528, algorithm 21.2)
% mean(Hq*X,2)=0 and svd(Hq*X)=1
HqS=normalizePoint(Hq*SOneClean,4);
c=mean(HqS,2);
HqSCentered=bsxfun(@minus,HqS,c);
% make sure the moments are 1
[U,S]=eig(HqSCentered*HqSCentered');
S=diag(1./sqrt(diag(S)));
K=S*U'; if det(K) < 0; K = -K; end
Ha=[ K -K*c; 0 0 0 1 ];
% define the quasi-affine inifinity that also centers the data
Hqa=Ha*Hq;

% Yalmip does not work under octave sorry :(
if exist('OCTAVE_VERSION','builtin')==5; return; end

% Chandraker IJCV 2009

% figure out the bounds on vi (equation 8)
PQuasi=multiTimes(P,inv(Hqa),1);
SQuasi=normalizePoint(Hqa*normalizePoint(SOneClean,-4),4);

% build the chirality inequalities
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
    if exitflag==1; if j==-1; u(i)=sol(i); else l(i)=sol(i); end; end
  end
end

[l,u]
%  A=zeros(3*nFrame-3,3); a=zeros(3*nFrame-3,1);
%  for i=2:nFrame
%    HInf=P(:,:,i)/eye(3,4);
%    A(3*i-5:3*i-3,:)=P(:,1:3,j)-HInf;
%    a(3*i-5:3*i-3,:)=P(:,4,j);
%  end
%  a\A
%  save('temp')
%  pause

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
yp=sdpvar(nFrame,2); t=sdpvar(nFrame,2); ypp=sdpvar(nFrame,2);

% rotated Lorentz cone constraint
FIni=set(cone([2*(fg(:,1)-fg(:,2));r-e],r+e)) + ...
  set(e(:)>=0) + set(r(:)>=0);

currBest=Inf; lowerBound=Inf;
eps=1e-4;
solverSetting = sdpsettings('solver','sdpa,csdp,sedumi,*','verbose',0, ...
  'cachesolvers',1);
for nItr=1:300
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

    % add relaxation constraints on f and g (ypp is for yp')
    for k=1:2
      F = F + relaxx13y(fg(:,k),greekHqat(:,:,k+2)*[v;1],...
        abcdBound(:,1,k+2),abcdBound(:,2,k+2),greekHqat(:,:,k)*[v;1],...
        abcdBound(:,1,k),abcdBound(:,2,k),yp(:,k),ypp(:,k),t(:,k));
    end
    % add the constraint on e
    F = F + concx83(e,greekHqat(1,:,4)*[v;1],abcdBound(1,1,4),...
      abcdBound(1,2,4));

    % solve the problem
    diagno = solvesdp(F,r,solverSetting);
    lowerBound(ind)=double(r);

    % make sure the best value is in the interval
    vBestTmp=min([max([double(v),l(:,ind)],[],2),u(:,ind)],[],2);

    % check for a few random values in the interval and check if they are
    % better than the value at the optimal v
    allV=[ vBestTmp, bsxfun(@plus,bsxfun(@times,rand(3,100),...
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

    if double(r) > currBest(ind); save('temp'); 'pprpprp', pause; end
  end
  % remove intervals for which the lower bound is higher than the current
  % best of another interval
  mini=min(currBest);
  badInterval=find(lowerBound>min(currBest));
  %    [l;u]
  % the following is necessary because of roundups
  if ~isempty(badInterval) && length(badInterval)~=length(currBest)
    l(:,badInterval)=[]; u(:,badInterval)=[]; vBest(:,badInterval)=[];
    currBest(badInterval)=[]; lowerBound(badInterval)=[];
  end
  lowerBound
  currBest
  [l;u]
  min(currBest)
  if size(l,2)==1 && norm(l-u,'fro')<eps; break; end
end
% figure out the best v
[disc,ind]=min(currBest);
vBest=vBest(:,ind);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function conc=concx83(z,x,xl,xu)
% compute the concave relaxation for x^(8/3)
conc=[ z <= nthroot(xl,3).^8+(x-xl)./(xu-xl).*(nthroot(xu,3).^8-...
  nthroot(xl,3).^8) ];
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function constr=relaxxyConv(z,x,xl,xu,y,yl,yu)
% compute the convex relaxation for x*y
constr=[ z >= xl.*y + yl.*x - xl.*yl; z >= xu.*y + yu.*x - xu.*yu ];
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function constr=relaxxyConc(z,x,xl,xu,y,yl,yu)
% compute the concave relaxation for x*y
constr=[ z <= xu.*y + yl.*x - xu.*yl; z <= xl.*y + yu.*x - xl.*yu ];
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function constr=relaxxy(z,x,xl,xu,y,yl,yu)
% precompute some date (because yalmip is slow ...)
xly=xl.*y; ylx=yl.*x; xuy=xu.*y; yux=yu.*x;
% compute the concave relaxation for x*y
constr=[ z <= xuy + ylx - xu.*yl; z <= xly + yux - xl.*yu ];
% compute the convex relaxation for x*y
constr=[ constr; z >= xly + ylx - xl.*yl; z >= xuy + yux - xu.*yu ];
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function obj=lineJoining(x,xl,xu)
% basically the line in equation A.10
obj=((nthroot(xu,3)-nthroot(xl,3))./(xu-xl)).*x+...
  (xu.*nthroot(xl,3)-xl.*nthroot(xu,3))./(xu-xl);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [obj, constr]=relaxx13yxposypos(x,xl,xu,y,yl,yu,yp)
% compute the convex relaxation for x^(1/3)*y when xl > 0 and yl > 0
% returns the under estimator and the two inequality constraints
obj=nthroot(xl,3).*yp+nthroot(xu,3).*(y-yp);
% precompute lam.*yl and lam.*yu (because yalmip is slow
lamyl=(yl./(xu-xl)).*x-xl.*yl./(xu-xl);
lamyu=(yu./(xu-xl)).*x-xl.*yu./(xu-xl);
constr=[ yl-lamyl<=yp<=yu-lamyu; lamyl<=(y-yp)<=lamyu ];
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function constr=relaxx13yxpos(z,x,xl,xu,y,yl,yu,yp,ypp,t)
% compute the convex relaxation for x^(1/3)*y when xl > 0
% convex relaxations for yl > 0 (like the paper)
ind=yl>0;
[obj, subConstr]=relaxx13yxposypos(x(ind),xl(ind),xu(ind),y(ind),yl(ind),...
  yu(ind),yp(ind));
constr=[ obj<=z(ind); subConstr ];
% convex relaxations for yl <= 0
ind=yl<0;
% add the constraints on t=x^(1/3) and z=ty
constr=constr+[ t(ind)>=lineJoining(x(ind),xl(ind),xu(ind)) ] + ...
  relaxxyConv(z(ind),t(ind),nthroot(xl(ind),3),...
  nthroot(xu(ind),3),y(ind),yl(ind),yu(ind));

% concave relaxations for yu < 0 (like the paper)
ind=yu<0;
[obj, subConstr]=relaxx13yxposypos(x(ind),xl(ind),xu(ind),-y(ind),-yu(ind),...
  -yl(ind),ypp(ind));
constr=constr+[ z(ind)<=-obj; subConstr ];
% concave relaxations for yu >= 0
ind=yu>=0;
% add the constraints on t=x^(1/3) and z=ty
constr=constr+[ t(ind)<=(x(ind)+2*xu(ind))./(3*nthroot(xu(ind),3).^2); ...
  t(ind)<=(x(ind)+2*xl(ind))./(3*nthroot(xl(ind),3).^2);] + ...
  relaxxyConc(z(ind),t(ind),nthroot(xl(ind),3),...
  nthroot(xu(ind),3),y(ind),yl(ind),yu(ind));


%  constr=[ t>=lineJoining(x,xl,xu) ] + ...
%    [ t<=(x+2*xu)./(3*nthroot(xu,3).^2); ...
%    t<=(x+2*xl)./(3*nthroot(xl,3).^2)] + ...
%    relaxxy(z,t,nthroot(xl,3),nthroot(xu,3),y,yl,yu);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [constr,indTot]=relaxx13y(z,x,xl,xu,y,yl,yu,yp,ypp,t)
% compute the convex relaxation for x^(1/3)*y
% case 1
indTot=[];
ind=(xl>0); indTot(end+1)=nnz(ind);
constr=relaxx13yxpos(z(ind),x(ind),xl(ind),xu(ind),y(ind),...
  yl(ind),yu(ind),yp(ind),ypp(ind),t(ind));

% case 1 bis
ind=(xu<0); indTot(end+1)=nnz(ind);
constr=constr+relaxx13yxpos(z(ind),-x(ind),-xu(ind),-xl(ind),-y(ind),...
  -yu(ind),-yl(ind),yp(ind),ypp(ind),t(ind));

% case 2
% t convex constraints
ind=(xl<=0) & (0<=xu) & (-xu/8>xl); indTot(end+1)=nnz(ind);
constr=[ constr; ...
  t(ind)>=(nthroot(xu(ind),-3).^2-2/3*...
  nthroot(xl(ind),3)./xu(ind)).*x(ind)+2/3*nthroot(xl(ind),3); ...
  t(ind)>=(x(ind)+2*xl(ind))./(3*nthroot(-xl(ind),3).^2) ];
% the following violates Cauchy conditions
%  constr=[ constr; ...
%    t(ind)>=(nthroot(xu(ind),3).^(-2)*4/3.*(x(ind)-xu(ind))+ ...
%    nthroot(xu(ind),3)); ...
%    t(ind)>=(x(ind)+2*xl(ind))./(3*nthroot(-xl(ind),3).^2) ];
% t convex constraints bis
ind=(xl<=0) & (0<=xu) & (-xu/8<=xl); indTot(end+1)=nnz(ind);
constr=[ constr; t(ind)>= lineJoining(x(ind),xl(ind),xu(ind)) ];
% t concave constraints
ind=(xl<=0) & (0<=xu) & (-xl/8<xu); indTot(end+1)=nnz(ind);
constr=[ constr; ...
  t(ind)<=(nthroot(xl(ind),-3).^2-2/3*nthroot(xu(ind),3)./xl(ind))...
  .*x(ind)+2/3*nthroot(xu(ind),3); ...
  t(ind)<=(x(ind)+2*xu(ind))./(3*nthroot(xu(ind),3).^2) ];
% the following violates Cauchy conditions
%  constr=[ constr; ...
%    t(ind)<=nthroot(xl(ind),3).^(-2)*4/3.*(x(ind)-xl(ind))+...
%    nthroot(xl(ind),3); ...
%    t(ind)<=(x(ind)+2*xu(ind))./(3*nthroot(xu(ind),3).^2) ];
% t concave constraints bis
ind=(xl<=0) & (0<=xu) & (-xl/8>=xu); indTot(end+1)=nnz(ind);
constr=[ constr; t(ind)<=lineJoining(x(ind),xl(ind),xu(ind)) ];

% add the constraints on z=ty
ind=(xl<=0) & (0<=xu); indTot(end+1)=nnz(ind);
constr=[ constr; relaxxy(z(ind),t(ind),nthroot(xl(ind),3),...
  nthroot(xu(ind),3),y(ind),yl(ind),yu(ind)) ];
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
