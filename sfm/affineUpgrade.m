function anim = affineUpgrade(anim)
% Perform an affine upgrade
%
% Given projection matrices and a 3D projective structure in anim,
% perform an affine upgrade.
% For now, only quasi-affine is performed
%
% USAGE
%   anim = affineUpgrade(anim, isCalibrated)
%
% INPUTS
%  anim       - Animation object with P and S filled
%
% OUTPUTS
%  anim      - Animation object after the quasi affine transformation
%              has been applied
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
    [sol,fval,exitflag] = linprog([0,0,0,0,-1]',-A,zeros(size(A,1),1),[],[],...
      [-1,-1,-1,-1,-Inf]', [1,1,1,1,Inf]');
    if exitflag==1 && sol(5)>0; break; end
  end
end
% define Hq
Hq=eye(4); Hq(4,:)=sol(1:4)'; Hq(1,1)=delta/Hq(4,4);

% define Ha
Sq=normalizePoint(Hq*normalizePoint(anim.S,-4),4);
Ha=eye(4); Ha(1:3,4)=-mean(Sq,2);
Hqa=Ha*Hq;

% apply Hqa
anim.P=multiTimes(anim.P,inv(Hqa),1);
anim.S=normalizePoint(Hqa*normalizePoint(anim.S,-4),4);

if true || exist('OCTAVE_VERSION','builtin')==5; return; end

% Chandraker IJCV 2009

r=sdpvar(1,1); e=sdpvar(1,1); f=sdpvar(nFrame,1); g=sdpvar(nFrame,1);
v=sdpvar(3,1);
% free variables for the convex/concave relaxations
ap=sdpvar(nFrame,1); lam=sdpvar(nFrame,1);

FIni=set(r*e >= sum((f-g).^2));

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
  vBest=zeros(s3,size(l,2)); currBest=zeros(1,size(l,2));
  for i=1:size(l,2)
    % add boundary conditions
    F=FIni+set(l(:,i)<=v)+set(v<=u(:,i));

    %  syms lam a1 a2 a3 a4 a5 a6 a7 a8 a9 a10 a11 a12 p1 p2 p3;
    %  A=[a1 a2 a3 a4; a5 a6 a7 a8; a9 a10 a11 a12];
    %  alpha = -(a4*p1 + a8*p2 + a12*p3 - a11 - a6 - a1);
    %  beta = (a2*a8 + a12*a3 - a11*a4  - a4*a6)*p1 + (a12*a7 - a1*a8 + a4*a5 - a11*a8)*p2 + (- a1*a12 - a12*a6 + a10*a8 + a4*a9)*p3 + a1*a11 + a1*a6 - a2*a5 + a11*a6 - a10*a7 - a3*a9;
    %  gamma = - ((a12*a2*a7 - a12*a3*a6 - a11*a2*a8 - a10*a4*a7 + a11*a4*a6 + a10*a3*a8)*p1 + (-a1*a12*a7 + a12*a3*a5 + a1*a11*a8 - a11*a4*a5 - a3*a8*a9 + a4*a7*a9)*p2 + (a1*a12*a6 - a12*a2*a5 - a1*a10*a8 + a10*a4*a5 + a2*a8*a9 - a4*a6*a9)*p3 + a11*a2*a5 - a1*a11*a6 + a1*a10*a7 - a10*a3*a5 - a2*a7*a9 + a3*a6*a9);
    %  simplify(lam^3 - alpha*lam^2 + beta*lam - gamma-det(lam*eye(3)-(A(:,1:3)-A(:,4)*[p1,p2,p3])))
    alpha=zeros(nFrame,j); beta=zeros(nFrame,j); gamma=zeros(nFrame,j);
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
    % greek will designate alpha,bet,gamma,delta
    greekHqat=zeros(nFrame,4,4);
    abcdBound=zeros(nFrame,2);
    for j=1:4
      switch j,
        case 1,
          greek=alpha;
      end
      greekHqat(:,:,j)=greek*Hqa';
      pos=greekHqat(:,:,j)>0; neg=greekHqat(:,:,j)<0;
      if j<4
        abcdBound(pos,:)=abcdBound(pos>0,:)+...
          greekHqat(pos,:,j)*[l(j,i),u(j,i)];
        abcdBound(neg,:)=abcdBound(neg,:)+...
          greekHqat(neg,:,j)*[u(j,i),l(j,i)];
      else
        abcdBound(pos,:)=abcdBound(pos,:)+[l(j,i),u(j,i)];
        abcdBound(neg,:)=abcdBound(neg,:)+[u(j,i),l(j,i)];
      end
    end
    % add first relaxation constraints
    ind1=abcdBound(:,1)>0 || abcdBound(:,2)<0;
    ind2=abcdBound(:,1)<0 && abcdBound(:,2)>0;
    % conv(c^(1/3)*a) <= f <= conc(c^(1/3)*a)
    F=F+set(f(ind1)>=(abcdBound(:,1,3).^(1/3)).*ap(ind1)+...
      (abcdBound(:,2,3).^(1/3))*(a(ind1)-ap(ind1)));
    F=F+set(f(ind1)>=(abcdBound(:,1,3).^(1/3)).*ap(ind1)+...
      (abcdBound(:,2,3).^(1/3))*(a(ind1)-ap(ind1)));
    % add second relaxation constraints
    % add third relaxation constraints
    % solve the problem
    diagno = solvesdp( F,r,sdpsettings('solver', ...
      'sdpa,csdp,sedumi,*','dualize',1,'debug',1));
    vBest(:,i)=double(v);
    % compute the current best
    a=alphaHqat*[v;1]; b=betaHqat*[v;1]; c=gammaHqat*[v;1];
    d=Hqa(:,4)'*[v;1];
    curBest(:,i)=sum((c^(1/3).*a-d^(1/3).*b)^2)/d^(8/3);
  end
end

end

function conc=concx83(z,x,xl,xu)
  % compute the concave relaxation for x^(8/3)
  conc=set(z <= xl.^(8/3)+(x-xl)./(xu-xl).*(xu.^(8/3)-xl.^(8/3)));
end

function conc=concxy(z,x,xl,xu,y,yl,yu)
  % compute the concave relaxation for x*y
  conc=set(z <= xu.*y + yl.*x - xu.*yl ) + ...
    set(z <= xl.*y + yu.*x - xl.*yu );
end

function conv=convxy(z,x,xl,xu,y,yl,yu)
  % compute the convex relaxation for x*y
  conc=set(z >= xl.*y + yl.*x - xl.*yl ) + ...
    set(z >= xu.*y + yu.*x - xu.*yu );
end

function conc=concx13y(z,x,xl,xu,y,yl,yu,ypp,lam,t)
  % compute the concave relaxation for x^(1/3)*y
  % case 1
  ind=(xl>0) || (xu<0);
  conc=set(z(ind)<=(xu(ind).^(1/3)-xl.^(1/3)).*ypp(ind) + ...
    xu(ind).^(1/3).*y(ind));
  conc=conc + set(lam(ind).*yu(ind) >= y(ind)+ypp(ind)) + ...
    set(y(ind)+ypp(ind)>=lam(ind).*yl(ind));
  conc=conc + set((1-lam(ind)).*yu(ind) >= -ypp(ind)) + ...
    set(-ypp(ind)>=(1-lam(ind)).*yl(ind));
  conc=conc + set(lam(ind)==(x(ind)-xl(ind))./(xu(ind)-xl(ind)));
  % case 2
  ind=(xl<=0) && (0<=xu) && (-xl/8<xu);
  conc=conc + set(t<=(xl.^(-2/3)-2/3*xu.^(1/3)./xl).*x+2/3*xu.^(1/3));
  conc=conc + set(t<=(x+2*xu)./(3*xu.^(2/3)));
  conc=conc + concxy(z(ind),t(ind),xl(ind).^(1/3),2/3*xu(ind).^(2/3),...
    y,yl,yu);
  % case 2 bis
  ind=(xl<=0) && (0<=xu) && (-xl/8>=xu);
  conc=conc + set(t<=((xu.^(1/3)-xl.^(1/3)).*x+...
    (xu.*xl.^(1/3)-xl.*xu.^(1/3)))./(xu-xl));
  conc=conc + concxy(z(ind),t(ind),xl(ind).^(1/3),xu(ind).^(1/3),...
    y,yl,yu);
end

function conv=convx13y(z,x,xl,xu,y,yl,yu,yp,lam)
  % compute the convex relaxation for x^(1/3)*y
  % case 1
  ind=(xl>0) || (xu<0);
  conv=set(z(ind)>=xl(ind).^(1/3).*yp(ind)+xu.^(1/3).* ...
    (y(ind)-yp(ind)));
  conv=conv + set((1-lam(ind)).*yl(ind) <= yp(ind)) + ...
    set(yp(ind)<=(1-lam(ind)).*yu(ind));
  conv=conv + set(lam(ind).*yl(ind) <= y(ind)-yp(ind)) + ...
    set(y(ind)-yp(ind)<=lam(ind).*yu(ind));
  conv=conv + set(lam(ind)==(x(ind)-xl(ind))./(xu(ind)-xl(ind)));
  % case 2
  ind=(xl<=0) && (0<=xu) && (-xu/8>xl);
  conc=conc + set(t>=(xu.^(-2/3)-2/3*xl.^(1/3)./xu).*x+2/3*xl.^(1/3));
  conc=conc + set(t>=(x+2*xl)./(3*(-xl).^(2/3)));
  conc=conc + convxy(z(ind),t(ind),zeros(len(ind),1),...
    max([xl(ind).^(1/3),xu(ind).^(1/3)],2),y,yl,yu);
  % case 2 bis
  ind=(xl<=0) && (0<=xu) && (-xu/8<=xl);
  conc=conc + set(t<=((xu.^(1/3)-xl.^(1/3)).*x+...
    (xu.*xl.^(1/3)-xl.*xu.^(1/3)))./(xu-xl));
  conc=conc + convxy(z(ind),t(ind),xl(ind).^(1/3),xu(ind).^(1/3),...
    y,yl,yu);
end
