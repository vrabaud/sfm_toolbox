addpath(genpath('~/SSFM'));
addpath(genpath('~/matlab'));
addpath('~/LSML/LSML');
rmpath(genpath('/users/u1/vrabaud/LSML/LSML/milf'));


nFrame=size(anim.W,3); nPoint=size(anim.W,2);
animGood2=animGood1;
S = animGood2.S; R = animGood2.R; t = animGood2.t; W = animGood2.W;
Pi=eye(2,3);


  temp = R(:,:,2:end)-R(:,:,1:end-1);
  res=[]; for i=1:length(temp); res(i)=norm(temp(:,:,i),'fro'); end
  plot(res);




% continue optimizing to remove outliers
lamt=1; lamS=1; epsR=0.07; epst=0.01; lamR=nPoint/epsR*0.001; dof=1;
ssfmError(Pi,W,S,R,t,lamt,lamR,lamS,epsR,epst);


% perform one averaging
animGood3 = animGood1;
S=animGood3.S; manifold=[];
X=reshape(S,[],nFrame);
%neigh=computeNeighbor(X,struct('k',5)); Nmat=neigh.Nmat;
temp = zeros(1,nFrame); temp(1:5)=1;
Nmat = toeplitz(temp);
[Chi3, Chi] = LSMLdenoise( X, Nmat, manifold, struct('show',1,'aveFlag',true) );
S=reshape(Chi,3,[],nFrame);
temp=squeeze(S(:,20,:));
plot3(temp(1,:),temp(2,:),temp(3,:),'.');
animGood4=animGood3; animGood4.S=S;

[res,animBad]=compareAnim(anim,8);
res1=compareAnim(animBad); [mean(res1),sum(res1)]
res2=compareAnim(animGood4); [mean(res2),sum(res2)]

% Continue optimizing but by uysing the d.o.f. constraint.
animGood4=animGood1;
S = animGood4.S; R = animGood4.R; t = animGood4.t; W = animGood4.W;
for j=2:nFrame; R(:,:,j)=rotationMatrix( R(:,:,j) ); end
%[ S R ] = ssfmChirality( S, R, t );

% slinky:
%epst=0.0; lamS=1; epsR=0.0; lamR=0.1; lamt=1;

% shark:
%epst=0.0; lamS=1; epsR=0.0; dof=2;

% coaster:
%epst=0.0; lamS=1; epsR=0.0;

% t=optimizet( Pi, W, S, R, t, lamt, epst );
% ssfmError(Pi,W,S,R,t,lamt,lamR,lamS,epsR,epst);
% [ R S ]=optimizeR(Pi,W,S,R,t,lamR,lamS,epsR);
% ssfmError(Pi,W,S,R,t,lamt,lamR,lamS,epsR,epst);
% S=optimizeS(Pi,W,S,R,t,lamS,dof, 'all');
% ssfmError(Pi,W,S,R,t,lamt,lamR,lamS,epsR,epst);

for i=1:1000
  'new denoising'
  t=ssfmOptimizet( isProj, W, S, R, t, lamt, epst );
  [ err errTot ] = computeSFMError( 'msfm', anim.isProj, 'W', W, 'S', ...
    S, 'R', R, 't', t, 'lamt', lamt, 'lamR', lamR, 'lamS', lamS, ...
    'epsR', epsR, 'epst', epst );
  errTot

  [ R S ]=optimizeR(Pi,W,S,R,t,lamR,lamS,epsR);
  [ err errTot ] = computeSFMError( 'msfm', anim.isProj, 'W', W, 'S', ...
    S, 'R', R, 't', t, 'lamt', lamt, 'lamR', lamR, 'lamS', lamS, ...
    'epsR', epsR, 'epst', epst );
  errTot

  temp=squeeze(S(:,20,:)); figure(2); clf;
  plot3(temp(1,:),temp(2,:),temp(3,:),'.');

  figure(3); clf;
  temp1 = R(:,:,2:end)-R(:,:,1:end-1); temp2 = t(:,2:end)-t(:,1:end-1);
  res=zeros(2,nFrame-1);
  for j=1:nFrame-1; res(:,j)=[norm(temp1(:,:,j),'fro') ...
      norm(temp2(:,j),'fro')]'; end
  plot(res(1,:),'b'); hold on; plot(res(2,:),'r');

  drawnow
  %S = optimizeS( Pi, W, S, R, t, lamS, dof, 'dof',SimMax<0.2  );
  S = optimizeS( Pi, W, S, R, t, lamS, dof, 'dof' );
  [ err errTot ] = computeSFMError( 'msfm', anim.isProj, 'W', W, 'S', ...
    S, 'R', R, 't', t, 'lamt', lamt, 'lamR', lamR, 'lamS', lamS, ...
    'epsR', epsR, 'epst', epst );
  errTot

  S=optimizeS(Pi,W,S,R,t,lamS,dof, 'outlier');
  [ err errTot ] = computeSFMError( 'msfm', anim.isProj, 'W', W, 'S', ...
    S, 'R', R, 't', t, 'lamt', lamt, 'lamR', lamR, 'lamS', lamS, ...
    'epsR', epsR, 'epst', epst );
  errTot

  temp=squeeze(S(:,20,:)); figure(2); clf;
  plot3(temp(1,:),temp(2,:),temp(3,:),'.');

  animGood5=animGood4; animGood5.S=S; animGood5.R=R; animGood5.t=t;
  animGood5.t(3,:)=1; animGood5.W=W;
  res2=compareAnim(anim,animGood5); [mean(res2(1,:)),mean(res2(2,:))]
  figure(1); plot(res2(2,:)); drawnow

  save temp
  diffS = S(:,:,2:end)-S(:,:,1:end-1);
  temp=diffS(:,:,2:end)-diffS(:,:,1:end-1);
  temp=squeeze(sum(sum(temp.^2,1),2)); figure(4); plot(temp);
end



playAnim(animGood5,struct('fps',10,'nCam',0));



animGood3=animGood2;
animGood3.S=S;
animGood3.R=R;
animGood3.t=t;
animGood3.t(3,:)=1;
animGood3.W=W;

res1=compareAnim(anim,animBad);
[mean(res1)]
res2=compareAnim(anim,animGood2);
mean(res2(3,:))
res2=compareAnim(anim,animGood5);
mean(res2(2,:))
plot(res2(2,:))


X=reshape(S,[],nFrame);
neigh=computeNeighbor(X,struct('k',10));
Nmat=neigh.Nmat;
[Chi3, Chi] = LSMLdenoise( X, Nmat, manifold, struct('show',1,'aveFlag',true) );
%[Chi3, Chi] = LSMLdenoise( X, Nmat, manifold, struct('show',1,'lambdaN',1) );
S=reshape(Chi,3,[],nFrame);

playAnim(animGood3,struct('fps',10,'nCam',1));

temp=squeeze(animGood3.S(:,20,:));
plot3(temp(1,:),temp(2,:),temp(3,:),'.');

temp=squeeze(animGood5.S(:,20,:));
plot3(temp(1,:),temp(2,:),temp(3,:),'.');


temp = R(:,:,2:end)-R(:,:,1:end-1); res=[];
for i=1:length(temp); res(i)=norm(temp(:,:,i),'fro'); end
plot(res)

temp = animGood1.R(:,:,2:end)-animGood1.R(:,:,1:end-1); res=[];
for i=1:length(temp); res(i)=norm(temp(:,:,i),'fro'); end
plot(res)



%% Visualize manifold
X=reshape(S,[], size(animGood5.S,3));
neigh=computeNeighbor(X,struct('k',10)); Nmat=neigh.Nmat;
manifold=fwLSMLoptimizeTH(X,Nmat,struct('lam',1e-1,'d',1,'rbfk',80,...
'show',0,'nitr',10,'nrestart',20));

i=20;
temp = LSMLcomputeH(X(:,i),manifold);
temp2 = squeeze(max(sum(temp.^2,1),[],2));
temp(:,:,temp2>mean(temp2)+0.6*std(temp2))=0;
visualizeManifold(S(:,:,i),1,struct('VF',reshape(temp,3,1,[])));

i=20
temp = LSMLcomputeH(X(:,:),manifold);
temp = 3*temp(3*i-2:3*i,:);
%temp2 = squeeze(max(sum(temp.^2,1),[],2));
%temp(:,:,temp2>mean(temp2)+0.6*std(temp2))=0;
visualizeManifold(squeeze(S(:,i,:)),1,struct('VF',reshape(temp,3,1,[])));




[Chi3, Chi] = LSMLdenoise( X, Nmat, manifold, struct('show',0,...
   'lambdaN',0.01,'nitr',10) );
visualizeManifold(squeeze(S(:,i,:)),1,struct('Chi3',squeeze(Chi3(:,:,3*i-2:3*i))'));




X=reshape(S,[],nFrame);
neigh=computeNeighbor(X,struct('k',10)); Nmat=neigh.Nmat;
[Chi3, Chi] = LSMLdenoise( X, Nmat, manifold, struct('show',1,'lambdaN',1) );
