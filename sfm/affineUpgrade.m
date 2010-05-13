function [ Hqa, vBest ] = affineUpgrade(anim)
% Perform an affine upgrade by recovering the plane at infinity
%
% Given projection matrices and a 3D projective structure in anim,
% perform an affine upgrade.
% For now, only quasi-affine is performed
%
% USAGE
%   anim = affineUpgrade(anim)
%
% INPUTS
%  anim          - Animation object with P and S filled
%
% OUTPUTS
%  Hqa           - quasi affine homography
%  vBest         - best v 
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

% remove points for which one w is negative
SOneClean=S; SOneClean(:,any(w<0,2)) = []; SOneClean(4,:)=1;

% find a solution to the chiral inequalities
for delta=[-1:2:1]
  % we want to minimize -d, with the parameters [v1,v2,v3,v4,d]
  A=[SOneClean';delta*C']; A(:,5)=-1;
  if exist('OCTAVE_VERSION','builtin')==5
    [sol,fval] = linprog([0,0,0,0,-1]',-A,zeros(size(A,1),1),zeros(0,5),...
      zeros(0,5),[-1,-1,-1,-1,-Inf]', [1,1,1,1,Inf]');
	if ~any(isna(sol)) && sol(5)>0; break; end
  else
    [sol,fval,exitflag] = linprog([0,0,0,0,-1]',-A,zeros(size(A,1),1),...
      [],[],[-1,-1,-1,-1,-Inf]', [1,1,1,1,Inf]', [], ...
      optimset('Display','off'));
    if exitflag==1 && sol(5)>0; break; end
  end
end
% define Hq
Hq=eye(4); Hq(4,:)=sol(1:4)'; Hq(1,1)=delta/Hq(4,4);

% define Ha
Sq=normalizePoint(Hq*normalizePoint(anim.S,-4),4);
Ha=eye(4); Ha(1:3,4)=-mean(Sq,2);
Hqa=Ha*Hq;

% Yalmip does not work under octave sorry :(
if exist('OCTAVE_VERSION','builtin')==5; return; end
return
% Chandraker IJCV 2009

% free variables for the convex/concave relaxations
r=sdpvar(1,1); e=sdpvar(1,1); fg=sdpvar(nFrame,2);
v=sdpvar(4,1);
% the constraint in re >= ... is replaced by (r+e)^2 >= ...
% as r and e are positive
% (looser but an SOCP cone)
FIni=set(cone(fg(:,1)-fg(:,2),r+e)) + set(v(4)==1);

l=repmat(-100,3,1); u=-l;
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
  vBest=zeros(4,size(l,2));
  currBest=zeros(1,size(l,2));
  lowerBound=zeros(1,size(l,2));
  for ind=1:size(l,2)
    % add boundary conditions
    F=FIni+set(l(:,ind)<=v(1:3)<=u(:,ind));

    %  syms lam a1 a2 a3 a4 a5 a6 a7 a8 a9 a10 a11 a12 p1 p2 p3;
    %  A=[a1 a2 a3 a4; a5 a6 a7 a8; a9 a10 a11 a12];
    %  alpha = -(a4*p1 + a8*p2 + a12*p3 - a11 - a6 - a1);
    %  beta = (a2*a8 + a12*a3 - a11*a4  - a4*a6)*p1 + (a12*a7 - a1*a8 + a4*a5 - a11*a8)*p2 + (- a1*a12 - a12*a6 + a10*a8 + a4*a9)*p3 + a1*a11 + a1*a6 - a2*a5 + a11*a6 - a10*a7 - a3*a9;
    %  gamma = - ((a12*a2*a7 - a12*a3*a6 - a11*a2*a8 - a10*a4*a7 + a11*a4*a6 + a10*a3*a8)*p1 + (-a1*a12*a7 + a12*a3*a5 + a1*a11*a8 - a11*a4*a5 - a3*a8*a9 + a4*a7*a9)*p2 + (a1*a12*a6 - a12*a2*a5 - a1*a10*a8 + a10*a4*a5 + a2*a8*a9 - a4*a6*a9)*p3 + a11*a2*a5 - a1*a11*a6 + a1*a10*a7 - a10*a3*a5 - a2*a7*a9 + a3*a6*a9);
    %  simplify(lam^3 - alpha*lam^2 + beta*lam - gamma-det(lam*eye(3)-(A(:,1:3)-A(:,4)*[p1,p2,p3])))
    alpha=zeros(nFrame,4); beta=zeros(nFrame,4); gamma=zeros(nFrame,4);
    for j=1:nFrame
      a1=P(1,1,j); a2=P(1,2,j); a3=P(1,3,j); a4=P(1,4,j);
      a5=P(2,1,j); a6=P(2,2,j); a7=P(2,3,j); a8=P(2,4,j);
      a9=P(3,1,j); a10=P(3,2,j); a11=P(3,3,j); a12=P(3,4,j);
      alpha(j,:)=-[ a4, a8, a12, - a11 - a6 - a1 ];
	  beta(j,:)=[a2*a8 + a12*a3 - a11*a4  - a4*a6, a12*a7 - a1*a8 + ...
	    a4*a5 - a11*a8, - a1*a12 - a12*a6 + a10*a8 + a4*a9, a1*a11 + ...
	    a1*a6 - a2*a5 + a11*a6 - a10*a7 - a3*a9 ];
      gamma(j,:)=-[ a12*a2*a7 - a12*a3*a6 - a11*a2*a8 - a10*a4*a7 + ...
        a11*a4*a6 + a10*a3*a8, -a1*a12*a7 + a12*a3*a5 + a1*a11*a8 - ...
        a11*a4*a5 - a3*a8*a9 + a4*a7*a9, a1*a12*a6 - a12*a2*a5 - ...
        a1*a10*a8 + a10*a4*a5 + a2*a8*a9 - a4*a6*a9, a11*a2*a5 - ...
        a1*a11*a6 + a1*a10*a7 - a10*a3*a5 - a2*a7*a9 + a3*a6*a9 ];
    end
    % compute the bounds on a,b,c,d
    % greek will designate alpha,bet,gamma
    greekHqat=zeros(nFrame,4,3);
    abcdBound=zeros(nFrame,2,4);
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
      % figure out the bounds of aj, bj, cj, d
      allProduct=bsxfun(@times,greekHqat(:,:,j),reshape(...
        [l(:,ind) u(:,ind); 1 1],1,4,2));
      abcdBound(:,:,j)=reshape(sum(sort(allProduct,3),2),nFrame,2);
    end
    abcdBound
    % add constraints
    lam=sdpvar(nFrame,2); yp=sdpvar(nFrame,2); ypp=sdpvar(nFrame,2);
	t=sdpvar(nFrame,2);
    for k=1:0
      F = F + convx13y(fg(:,k),greekHqat(:,:,k+2)*v,abcdBound(:,1,k+2),...
        abcdBound(:,2,k+2),greekHqat(:,:,k+1)*v,abcdBound(:,1,k+1),...
        abcdBound(:,2,k+1),yp(:,k),lam(:,k),t(:,k));
      F = F + concx13y(fg(:,k),greekHqat(:,:,k+2)*v,abcdBound(:,1,k+2),...
        abcdBound(:,2,k+2),greekHqat(:,:,k+1)*v,abcdBound(:,1,k+1),...
        abcdBound(:,2,k+1),ypp(:,k),lam(:,k),t(:,k));
    end
    F = F + concx83(e,greekHqat(:,:,4)*v,abcdBound(:,1,4),abcdBound(:,2,4));
    for j=1:length(F)
    F{j};
    end
    % solve the problem
    diagno = solvesdp( F,r,sdpsettings('solver', ...
      'sdpa,csdp,sedumi,*','dualize',true,'verbose',1));
    vBest(:,ind)=double(v);
    % compute the current best
    a=greekHqat(:,:,1)*vBest(:,ind); b=greekHqat(:,:,2)*vBest(:,ind);
    c=greekHqat(:,:,3)*vBest(:,ind); d=greekHqat(1,:,4)*vBest(:,ind);
    curBest(ind)=sum((nthroot(c,3).*a-nthroot(d,3).*b).^2)./nthroot(d,3)^8;
    lowerBound(ind)=double(r);
  end
  % remove intervals for which the lower bound is higher than the current
  % best of another intervals
  badInterval=find(lowerBound>min(curBest));
  [l;u]
  lowerBound
  curBest
  l(:,badInterval)=[]; u(:,badInterval)=[]; vBest(:,badInterval)=[];
  [l;u]
  'toto'
  
end
% figure out the best v
vBest=v(:,1);
end

function conc=concx83(z,x,xl,xu)
  % compute the concave relaxation for x^(8/3)
  conc=set(z <= nthroot(xl,3).^8+(x-xl)./(xu-xl).*(nthroot(xu,3).^8-...
    nthroot(xl,3).^8));
end

function conc=concxy(z,x,xl,xu,y,yl,yu)
  % compute the concave relaxation for x*y
  conc=set(z <= xu.*y + yl.*x - xu.*yl ) + ...
    set(z <= xl.*y + yu.*x - xl.*yu );
end

function conv=convxy(z,x,xl,xu,y,yl,yu)
  % compute the convex relaxation for x*y
  conv=set(z >= xl.*y + yl.*x - xl.*yl ) + ...
    set(z >= xu.*y + yu.*x - xu.*yu );
end

function relax=relaxx13y(z,x,xl,xu,y,yl,yu,yp,lam,isConv)
  % compute the convex relaxation for x^(1/3)*y when xl > 0
  % case 1
  if isConv
    relax=set(z>=nthroot(xl,3).*yp+nthroot(xu,3).*(y-yp));
  else
    relax=set(z<=nthroot(xl,3).*yp+nthroot(xu,3).*(y-yp));
  end
  relax=relax + set((1-lam).*yl<=yp<=(1-lam).*yu);
  relax=relax + set(lam.*yl<=y-yp<=lam.*yu);
  relax=relax + set(lam==(x-xl)./(xu-xl));
end

function conc=concx13y(z,x,xl,xu,y,yl,yu,ypp,lam,t)
  % compute the concave relaxation for x^(1/3)*y
  % case 1
  ind=(xl>0);
  conc=relaxx13y(z(ind),x(ind),xl(ind),xu(ind),y(ind),yl(ind),yu(ind),...
    ypp(ind),lam(ind),false);
  ind=(xu<0);
  conc=conc+relaxx13y(z(ind),x(ind),-xu(ind),-xl(ind),y(ind),yl(ind),...
    yu(ind),ypp(ind),lam(ind),false);
  % case 2
  ind=(xl<=0) & (0<=xu) & (-xl/8<xu);
  conc=conc + set(t(ind)<=(nthroot(xl(ind),-3).^2-2/3*nthroot(xu(ind),3)./xl(ind))...
    .*x(ind)+2/3*nthroot(xu(ind),3));
  conc=conc + set(t(ind)<=(x(ind)+2*xu(ind))./(3*nthroot(xu(ind),3).^2));
  conc=conc + concxy(z(ind),t(ind),nthroot(xl(ind),3),2/3*nthroot(xu(ind),3).^2,...
    y(ind),yl(ind),yu(ind));
  % case 2 bis
  ind=(xl<=0) & (0<=xu) & (-xl/8>=xu);
  conc=conc + set(t(ind)<=((nthroot(xu(ind),3)-nthroot(xl(ind),3)).*x(ind)+...
    (xu(ind).*nthroot(xl(ind),3)-xl(ind).*nthroot(xu(ind),3)))./(xu(ind)-xl(ind)));
  conc=conc + concxy(z(ind),t(ind),nthroot(xl(ind),3),nthroot(xu(ind),3),...
    y(ind),yl(ind),yu(ind));
  % add the constraint on t
  ind=(xl<=0) & (0<=xu);
  tBound=[nthroot(xl(ind),3),2/3*nthroot(xu(ind),3).^2, nthroot(xu(ind),3)];
  conc=conc+concxy(z(ind),t(ind),min(tBound,[],2),max(tBound,[],2),...
    y(ind),yl(ind),yu(ind));
end

function conv=convx13y(z,x,xl,xu,y,yl,yu,yp,lam,t)
  % compute the convex relaxation for x^(1/3)*y
  % case 1
  ind=(xl>0);
  conv=relaxx13y(z(ind),x(ind),xl(ind),xu(ind),y(ind),yl(ind),yu(ind),...
    yp(ind),lam(ind),true);
  ind=(xu<0);
  conv=conv+relaxx13y(z(ind),x(ind),-xu(ind),-xl(ind),y(ind),yl(ind),...
    yu(ind),yp(ind),lam(ind),true);
  % case 2
  ind=(xl<=0) & (0<=xu) & (-xu/8>xl);
  conv=conv + set(t(ind)>=(nthroot(xu(ind),-3).^2-2/3*...
    nthroot(xl(ind),3)./xu(ind)).*x(ind)+2/3*nthroot(xl(ind),3));
  conv=conv + set(t(ind)>=(x(ind)+2*xl(ind))./(3*nthroot(-xl(ind),3).^2));
  conv=conv + convxy(z(ind),t(ind),zeros(nnz(ind),1),...
    max([nthroot(xl(ind),3),nthroot(xu(ind),3)],[],2),y(ind),yl(ind),yu(ind));
  % case 2 bis
  ind=(xl<=0) & (0<=xu) & (-xu/8<=xl);
  conv=conv + set(t(ind)<=((nthroot(xu(ind),3)-nthroot(xl(ind),3)).*x(ind)+...
    (xu(ind).*nthroot(xl(ind),3)-xl(ind).*nthroot(xu(ind),3)))./(xu(ind)-xl(ind)));
  conv=conv + convxy(z(ind),t(ind),nthroot(xl(ind),3),nthroot(xu(ind),3),...
    y(ind),yl(ind),yu(ind));
  % add the constraint on t
  ind=(xl<=0) & (0<=xu);
  tBound=[nthroot(xl(ind),3),2/3*nthroot(xu(ind),3).^2, nthroot(xu(ind),3)];
  conv=conv+convxy(z(ind),t(ind),min(tBound,[],2),max(tBound,[],2),...
    y(ind),yl(ind),yu(ind));
end
