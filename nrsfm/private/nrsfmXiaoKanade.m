function anim = nrsfmXiaoKanade( W, nBasis )
% Compute orthographic non-rigid structure from motion using Xiao-Kanade
%
% Just check the Xiao's CVPR04 paper
%
% USAGE
%  anim = nrsfmXiaoKanade( W, nBasis )
%
% INPUTS
%  W             - [ 2 x nPoint x nFrame ] set of 2D points
%  nBasis        - {Computed } number of bases to use
%
% OUTPUTS
%  anim          - Animation object (help Animation for details)
%
% EXAMPLE
%
% See also COMPUTESMFROMW
%
% Vincent's Structure From Motion Toolbox      Version 3.0
% Copyright (C) 2009 Vincent Rabaud.  [vrabaud-at-cs.ucsd.edu]
% Please email me if you find bugs, or have suggestions or questions!
% Licensed under the GPL [see external/gpl.txt]

nFrame = size( W, 3 ); nPoint = size( W, 2 );
animBest=Animation; animBest.W = W;

W = reshape( permute( W, [ 1 3 2 ] ), [], nPoint );
animBest.t = reshape( mean(W,2), 2, nFrame ); animBest.t(3,:) = 0;
W = W - repmat( mean(W,2), [ 1 nPoint ] );

if nargin<2 || isempty(nBasis)
  [ U S V ] = svd(W,'econ'); S = diag(S); sumVal = S(1); Kd = 1;
  while ( Kd<=length(S) ) && ( sumVal < 0.995*sum(abs(S)) )
    if S(Kd)<0; break; end
    Kd = Kd + 1; sumVal = sumVal + S(Kd);
  end
  W = U(:,1:Kd)*diag(S(1:Kd))*V(:,1:Kd)';
else
  Kd = 3*nBasis;
  % Check if Kd is not smaller (degenerate deformations)
  [ U S V ] = svd(W,'econ');
  S=cumsum(diag(S));
  ind = find(S>=0.995*S(end));
  Kd=min(Kd,ind(1));
end

%%%%%%%%%%%%%%%%%%%%%% Determine K3 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
K3 = floor(Kd/3);
% Basis constraints
WTildeOri = [ W(1:2:end,:) W(2:2:end,:) ];

while 1
  % Perform random search for 5 seconds
  t0 = clock; mini = Inf;
  while etime( clock, t0 )<5
    indWTmp = randSample(nFrame,K3)';
    tmp = cond( WTildeOri( indWTmp, : ) );
    if tmp < mini; mini = tmp; indW = indWTmp; end
  end
  
  if mini>1e2 && K3>0; K3=K3-1; else break; end
end

%%%%%%%%%%%%%%%%%%%%%% Deal with the K3 basis %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Compute initial variables
[ MTilde S disc ] = svd( W );
MTilde = MTilde( :, 1:Kd )*S(1:Kd,1:Kd); gj3=zeros(Kd,3,K3);

% Rotation constraints
objRotation = 0; Qi = sdpvar( Kd, Kd );
F = set( Qi >= 0 );
for i = 1 : nFrame
  objRotation = objRotation + abs( MTilde(2*i-1,:)*Qi*MTilde(2*i-1,:)' - ...
    MTilde(2*i,:)*Qi*MTilde(2*i,:)' ) + ...
    abs( MTilde(2*i,:)*Qi*MTilde(2*i-1,:)' );
end

% Solve for Qk using all the constraints
kk = 1;
for i = indW
  nc = 2*nFrame + 1;
  objBasis = 0;
  disp('Computing the basis constraints');
  for m = indW
    if m==i; interN = m; else interN = 1 : nFrame; end
    for n = interN % could be sampled
      if m==n && m==i
        objBasis = objBasis +abs(MTilde(2*m-1,:)*Qi*MTilde(2*n-1,:)'-1) ...
          + abs( MTilde(2*m,:)*Qi*MTilde(2*n,:)' - 1);
      else
        objBasis = objBasis +abs(MTilde(2*m-1,:)*Qi*MTilde(2*n-1,:)') ...
          + abs( MTilde(2*m,:)*Qi*MTilde(2*n,:)' );
      end
      
      objBasis = objBasis + abs( MTilde(2*m-1,:)*Qi*MTilde(2*n,:)' ) + ...
        abs( MTilde(2*m,:)*Qi*MTilde(2*n-1,:)' );
    end
  end
  
  % Sample the constraints
  disp('Solving the SDP');
  diagno = solvesdp( F,objRotation+objBasis,sdpsettings('solver',...
    'sdpa,csdp,sedumi,*','dualize',1,'debug',1));
  
  QiTmp = double(Qi);
  
  [ U S disc ] = svd( QiTmp ); gj3(:,:,kk) = U(:,1:3)*sqrt(S(1:3,1:3));
  kk = kk + 1;
end

%  [U S V]=svd(constr,'econ');
%  rank(constr)
%  size(null(constr(1:2*nFrame,:)))
%  size(constr)
%  tmp=diag(S)';
%  tmp(end-3:end)
%  size(null(constr),2)
%  return
%%%%%%%%%%%%%%%%%%%%%% Get K2 and K1 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Compute K2
K2=floor((Kd-3*K3)/2);
while K2>0
  if (K2^2+K2)/2>size(null(constr),2); K2=K2-1; else break; end
end

%% The explanations given by Xiao04 seem wrong to determine K1 and K2
%% So, let's do it the lazy way: let's try all possible combinations of K1
%% and K2
errBest = Inf;
fprintf('\nSo far, Kd=%d, K3=%d\n', Kd, K3 );
for K2 = floor((Kd-3*K3)/2) : -1 : 0
  K1 = Kd - 3*K3 - 2*K2;
  K = K3+K2+K1;
  
  fprintf('\nTrying, Kd=%d, K1=%d, K2=%d, K3=%d\n', Kd, K1, K2, K3 );
  
  %%%%%%%%%%%%%%%%%%%%%% Get the sets of rotations %%%%%%%%%%%%%%%%%%%%%%%%
  G = reshape( gj3, Kd, [] );
  M = MTilde*G; M = M/max(abs(M(:)))*K3; % Just for numerical stability
  R = zeros( 2*nFrame, 3*K3 ); l = zeros( K, nFrame ); thres = 0.001;
  rotIsBad = zeros( K3, nFrame );
  for i = 1 : nFrame
    for k = 1 : K3
      tmp = M(2*i-1:2*i,3*k-2:3*k);
      l(k,i) = norm(tmp);
      
      if abs( l(k,i) )>thres; R(2*i-1:2*i,3*k-2:3*k) = tmp/l(k,i);
      else rotIsBad(k,i)=1;
      end
    end
  end
  
  % Rectify R (not really mentioned by Xiao)
  for k=2:K3
    for i = 2:nFrame
      goodFrame = find( sum( rotIsBad([1 k], 1:i), 1 )==0);
      goodFrame = sort( [ 2*goodFrame 2*goodFrame-1 ] );
      if length(goodFrame)<4; continue; end
      R1 = R( goodFrame, 1:3 );
      R2 = R( goodFrame, 3*k-2:3*k );
      R3 = R2; R3(end-1:end,:) = -R3(end-1:end,:);
      
      if norm( R2*(R2\R1)-R1, 'fro' )>norm( R3*(R3\R1)-R1, 'fro' )
        R( goodFrame(end-1:end), 3*k-2:3*k ) = ...
          -R( goodFrame(end-1:end), 3*k-2:3*k );
      end
    end
    
    goodFrame = find( sum( rotIsBad([1 k], 1:nFrame), 1 )==0);
    goodFrame = sort( [ 2*goodFrame 2*goodFrame-1 ] );
    R1 = R( goodFrame, 1:3 );
    R2 = R( goodFrame, 3*k-2:3*k );
    tmp=R2\R1;
    
    R( :, 3*k-2:3*k ) = R( :, 3*k-2:3*k )*tmp;
    M( :, 3*k-2:3*k ) = M( :, 3*k-2:3*k )*tmp;
  end
  
  % Get all the rotations and coefficients in K3 (if non-degenerate)
  RTot = zeros(3,3,nFrame);
  for i = 1 : nFrame
    % Get the average rotation matrix
    tmp=find(~rotIsBad(:,i)');
    if isempty(tmp); continue; end
    RTmp=zeros(3,3,length(tmp));
    
    qTot=cell(1,2);
    for k=1:2
      kk=1;
      for j=tmp
        RTmp(:,:,kk)=rotationMatrix( R(2*i-1:2*i,3*j-2:3*j)*(2*k-3) );
        kk=kk+1;
      end
      
      q=quaternion(RTmp);
      for j=2:size(q,2)
        if norm(-q(:,j)-q(:,1))<norm(q(:,j)-q(:,1)); q(:,j)=-q(:,j); end
      end
      qTot{k}=mean(q,2);
    end
    
    q=qTot{1};
    if i>=2
      q0 = quaternion( RTot(:,:,i-1) );
      if min( norm(qTot{2}-q0), norm(-qTot{2}-q0) )<...
          min( norm(qTot{1}-q0), norm(-qTot{1}-q0) )
        q=qTot{2};
      end
    end
    
    RTot(:,:,i) = quaternion(q);
    
    % Recover the coefficients
    l(1:K3,i)=(vect(RTot(1:2,:,i),'v')\reshape(M(2*i-1:2*i,1:3*K3),6,K3))';
    
    % Compute MTilde
    MTilde(2*i-1:2*i,1:3*K3)=kron(l(1:K3,i)',RTot(1:2,:,i));
  end
  for k = 1 : K3
    gj3(:,:,k) = 0; gj3(3*k-2:3*k,:,k) = eye(3);
  end
  R = RTot;
  rotIsUnknown = sum(abs(l(1:K3,:))<thres,1)==K3;
  if sum(~rotIsUnknown)==0
    [disc indMax]=max(sum(l(1:K3,:).^2,1)); rotIsUnknown(indMax)=1;
  end
  goodFrame = find( ~rotIsUnknown );
  tmp=sort([2*goodFrame,2*goodFrame-1]);
  
  % Get the MTilde where l(1:K3,i) is close to 0
  M=MTilde(tmp,:); B=M\W(tmp,:); M2=W/B;
  for i=find(rotIsUnknown)
    MTilde(2*i-1:2*i,:)=M2(2*i-1:2*i,:);
  end
  
  %   find(rotIsUnknown)
  %   plot(l')
  %     B=MTilde\W;
  %   norm(MTilde*B-W,'fro')
  %   pause
  
  % The first 3*K3 columns of MTilde and MHat are now identical
  % We now need to deal with the last columns
  
  %%%%%%%%%%%%%%%%%%%%%% Find gj and rj %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % We here perform a full optimization and not the hacky alternate
  % optimization from the paper
  MNull = zeros( 2*K3, Kd );
  for i = 1 : K3
    MNull(2*i-1:2*i,:) = MTilde( 2*indW(i)-1:2*indW(i), : );
  end
  [ U S V ] = svd( MNull, 'econ' );
  [ V disc ] =rq(V);
  gj2 = reshape( V(:,end-(2*K2+K1)+1:end-K1), [], 2, K2);rj2=zeros(3,2,K2);
  gj1 = V(:,end-K1+1:end); rj1=zeros(3,K1);
  
  % Perform iterations for K2 bases
  opt = optimset( 'Display', 'on', 'GradObj', 'off', 'Hessian', 'off',...
    'LargeScale', 'off' );
  for j=1:K2
    % Get the best rj2
    ATmp=zeros(4*length(goodFrame),3);
    for kk=1:2
      n = 1;
      for m=goodFrame
        for k=1:2
          ATmp(n,:)=gj2(:,k,j)'*MTilde(2*m-1,:)'*R(2,:,m) - ...
            gj2(:,k,j)'*MTilde(2*m,:)'*R(1,:,m);
          n = n + 1;
        end
      end
      
      % Force the rj2 to be independent
      zeroStart=2+kk;
      [ U S V ] = svd( ATmp(:,1:zeroStart-1), 'econ' );
      rj2(1:zeroStart-1,kk,j) = V(:,end)/V(end,end);
      rj2(zeroStart:end,kk,j) = 0;
    end
    
    % Perform the full optimization
    onePos=3*K3+1+2*(j-1);
    x=[gj2(1:onePos-1,1,j); gj2(1:onePos,2,j); rj2(1:2,1,j); rj2(1:2,2,j)];
    MTildeTmp= MTilde(:,1:onePos+1); MNullTmp=MNull(:,1:onePos+1);
    x = fminunc( @(X)( optimizeGj2( X,goodFrame, MTildeTmp, ...
      R, MNullTmp ) ), x,opt);
    
    % Save the best results
    rj2(:,:,j)=[x(end-3:end-2) x(end-1:end); 0 1 ];
    gj2(:,:,j)=0;
    gj2(1:onePos,1,j) = [ x(1:onePos-1); 1 ];
    gj2(1:onePos+1,2,j) = [ x(1:onePos); 1 ];
  end
  
  % Perform iterations for K1 bases
  for j=1:K1
    % Get the best rj1
    ATmp=zeros(length(goodFrame),3);
    for m=goodFrame
      ATmp(m,:) = gj1(:,j)'*( MTilde( 2*m-1 , : )'*R(2,:,m) - ...
        MTilde( 2*m , : )'*R(1,:,m) );
    end
    
    [ U S V ] = svd( ATmp, 'econ' );
    rj1(:,j) = V(:,end);
    
    % Perform the full optimization
    onePos=3*K3+2*K2+j;
    x=[gj1(1:onePos-1,j); rj1(:,j) ];
    MTildeTmp= MTilde(:,1:onePos); MNullTmp=MNull(:,1:onePos);
    x = fminunc( @(X)( optimizeGj1( X,goodFrame, MTildeTmp, ...
      R, MNullTmp ) ), x,opt);
    
    % Save the best results
    rj1(:,j)=x(end-2:end);
    gj1(:,j)=0; gj1(1:onePos,j) = [ x(1:onePos-1); 1 ];
  end
  
  % % quality
  % err=zeros(1,length(goodFrame));
  % for j = 1:K2
  %   for m=goodFrame
  %     err(m)=err(m)+abs(gj2(:,1,j)'*( MTilde( 2*m-1 , : )'*R(2,:,m) - ...
  %       MTilde( 2*m , : )'*R(1,:,m) )*rj2(:,1,j)) ...
  %        + abs(gj2(:,2,j)'*( MTilde( 2*m-1 , : )'*R(2,:,m) - ...
  %       MTilde( 2*m , : )'*R(1,:,m) )*rj2(:,2,j));
  %   end
  % end
  %       norm(err)
  % %       norm(MNull*gj1(:,j))
  %       plot(err)
  %
  % pause
  
  %%%%%%%%%%%%%%%%%%%%%% Recover the coefficients for K2 K1 %%%%%%%%%%%%%%%
  % Rectify the right columns of MTilde
  for m=1:nFrame
    for i = 1 : K2
      MTilde(2*m-1:2*m,3*K3+2*i-1:3*K3+2*i)=MTilde(2*m-1:2*m,:)*gj2(:,:,i);
    end
    for i = 1 : K1
      MTilde(2*m-1:2*m,3*K3+2*K2+i)=MTilde(2*m-1:2*m,:)*gj1(:,i);
    end
  end
  %   l(K3+i,m)*R(1:2,:,m)*rj2(:,:,i)=
  %                   MTilde(2*m-1:2*m,3*K3+2*i-1:3*K3+2*i)*RAmbiguityi
  %   l(K3+i,m)*R(1:2,:,m)*rj1(:,i)=
  %                   MTilde(2*m-1:2*m,3*K3+2*K2+i)
  for m=goodFrame
    for i = 1 : K2
      tmp1=R(1:2,:,m)*rj2(:,:,i); tmp1=tmp1*tmp1';
      tmp2=MTilde(2*m-1:2*m,3*K3+2*i-1:3*K3+2*i); tmp2=tmp2*tmp2';
      l(K3+i,m) = sqrt( tmp1(:)\tmp2(:) );
    end
    for i = 1 : K1
      tmp1=R(1:2,:,m)*rj1(:,i);
      tmp2=MTilde(2*m-1:2*m,3*K3+2*K2+i);
      l(K3+K2+i,m) = tmp1\tmp2;
    end
  end
  
  %   B=MTilde\W;
  %   norm(MTilde*B-W,'fro')
  %   plot(l')
  %   pause
  
  
  % Change the signs of the coefficients
  for i = 1 : K2
    A=zeros(2*nFrame,2);
    A(1:2,:)=l(K3+i,goodFrame(1))*R(1:2,:,goodFrame(1))*rj2(:,:,i);
    
    for j=2:length(goodFrame)
      m=goodFrame(j); A(2*j-1:2*j,:)=l(K3+i,m)*R(1:2,:,m)*rj2(:,:,i);
      tmp=goodFrame(1:j); tmp=sort( [ 2*tmp-1, 2*tmp ] );
      MTildeLoc=MTilde(tmp,3*K3+2*i-1:3*K3+2*i);
      
      RAmb1=MTildeLoc\A(1:2*j,:);
      err1=norm(MTildeLoc*RAmb1-A(1:2*j,:),'fro');
      
      A(2*j-1:2*j,:)=-A(2*j-1:2*j,:);
      RAmb2=MTildeLoc\A(1:2*j,:);
      err2=norm(MTildeLoc*RAmb2-A(1:2*j,:),'fro');
      
      if err2<err1; l(K3+i,m)=-l(K3+i,m);
      else A(2*j-1:2*j,:)=-A(2*j-1:2*j,:);
      end
    end
  end
  
  B=MTilde\W;
  %   norm(MTilde*B-W,'fro')
  %   plot(l')
  %   pause
  
  %%%%%%%%%%%%%%%%%%%%%% Recover the shape basis %%%%%%%%%%%%%%%%%%%%%%%%%%
  M=zeros(2*length(goodFrame),3*K);
  for i=1:length(goodFrame)
    M(2*i-1:2*i,:)=kron(l(:,goodFrame(i))',R(1:2,:,goodFrame(i)));
  end
  
  tmp=sort([2*goodFrame,2*goodFrame-1]);
  B=M\W(tmp,:);
  %   norm(M*B-W(tmp,:),'fro')
  %   tmp=sum((M*B-W(tmp,:)).^2,2);
  %   plot( tmp(1:2:end)+tmp(2:2:end) )
  %   'coiin'
  %   pause
  
  % Recover the full basis
  SBasis = zeros( 3, nPoint, K );
  for i = 1 : K; SBasis(:,:,i) = B(3*i-2:3*i,:); end
  
  %%%%%%%%%%%%%%%%%%%%%% Recover the missing coefficients %%%%%%%%%%%%%%%%%
  opt = optimset( 'Display', 'off', 'GradObj', 'off', 'Hessian', 'off',...
    'LargeScale', 'off' );
  
  while 1
    for m=1:nFrame
      if rotIsUnknown(m)
        if m==1 || rotIsUnknown(m-1)
          if m==nFrame; continue
          else if ~rotIsUnknown(m+1); iIni=m+1; else continue; end
          end
        else iIni=m-1;
        end
        
        % Optimize over the shape coefficients/rotation
        coeff=[ l(:,iIni); quaternion( R(:,:,iIni) ) ];
        Wm=W(2*m-1:2*m,:);
        
        [ c res ] = fminunc( @(X)( optimRl( X, K,B,Wm ) ), coeff, opt );
        
        l(:,m)=c(1:K); R(:,:,m)=quaternion(c(K+1:end));
        rotIsUnknown(m)=0;
      else
        if m>=2
          %           plot(quaternion(R)')
          %           pause
          q0 = quaternion( R(:,:,m-1) );
          q1 = quaternion( R(:,:,m) );
          q2 = quaternion( rotationMatrix( -R(1:2,:,m) ) );
          if min( norm(q2-q0), norm(-q2-q0) )<min(norm(q1-q0),norm(-q1-q0))
            l(:,m)=-l(:,m); R(:,:,m)=rotationMatrix(-R(1:2,:,m));
          end
        end
      end
    end
    if ~rotIsUnknown(1); break; end
  end
  
  %%%%%%%%%%%%%%%%%%%%%% Finalize the anim object %%%%%%%%%%%%%%%%%%%%%%%%%
  anim=Animation('SBasis', SBasis, 'l', l, 'R', R, 'isProj', false, ...
    't', animBest.t, 'W', animBest.W );
  err = anim.computeError( anim ); err=err(1);
  fprintf( 'Reprojection error: %f \n', err );
  
  if err<errBest; errBest=err; animBest=anim; end
end
anim=animBest;

fprintf( 'Best reprojection error: %f \n', errBest );
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function res = optimRl(X,K,B,Wm)
lm=X(1:K);
Rm=quaternion(X(K+1:end)/norm(X(K+1:end)));
res = norm(kron(lm',Rm(1:2,:))*B-Wm,'fro')^2;
end

function val = optimizeGj2( x, goodFrame, MTilde, R, MNull )
% gjk' MTilde_{2m-1}' R_m2 rjkk - gjk' MTilde_2m' R_m1 rjkk = 0,
% k and kk in {1 2}
% rj can be multiplied by any 2x2
rj2=[x(end-3:end-2) x(end-1:end); 0 1 ]; rj2(:,1)=rj2(:,1)/norm(rj2(:,1));
n=(length(x)-4-1)/2;
gj2=[ [ x(1:n); 1; 0 ] [ x(n+1:2*n+1); 1 ] ];
n=n+2;

n=1;
val=zeros(1,length(goodFrame));
A=gj2'*MTilde';
for m=goodFrame
  val(n)=norm( ( A(:,2*m-1,:)*R(2,:,m) - A(:,2*m)*R(1,:,m) )*rj2, 'fro' );
  n=n+1;
end
val=norm([ val vect(MNull*gj2,'v')' ]);
end

function val = optimizeGj1( x, goodFrame, MTilde, R, MNull )
% gj' ( MTilde_{2m-1} R_m2 - MTilde_2m R_m1 ) rj = 0
gj1=[ x(1:end-3); 1 ];
rj1=x(end-2:end); rj1=rj1/norm(rj1);

n=1;
A=gj1'*MTilde';
B=zeros(length(goodFrame),3);
for m=goodFrame
  B(n,:)=A(:,2*m-1)*R(2,:,m) - A(:,2*m)*R(1,:,m);
  n=n+1;
end
val=norm([ (B*rj1); vect(MNull*gj1,'v') ]);

end
