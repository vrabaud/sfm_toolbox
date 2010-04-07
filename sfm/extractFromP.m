function [K,R,t]=extractFromP(P,isProj,isKIdentity)
% Extract information from a projection matrix
%
% USAGE
%  [K,R,t]=extractFromP(P,isProj)
%
% INPUTS
%  P       - P, projection matrix
%  isProj  - flag indicating if the camera ia projective one
%
% OUTPUTS
%  K      - intrinsic parameter matrix
%  R      - rotation matrix
%  t      - translation vector
%
% EXAMPLE
%
% See also
%
% Vincent's Structure From Motion Toolbox      Version 3.0
% Copyright (C) 2009 Vincent Rabaud.  [vrabaud-at-cs.ucsd.edu]
% Please email me if you find bugs, or have suggestions or questions!
% Licensed under the GPL [see external/gpl.txt]

if nargin<2; isKIdentity=false; end

if isKIdentity
  return
end
if isProj
  if nargout==1
    % Use Kruppa equations to extract f
    % Reference: HZ2, p472, Example 19.8
    F=convertPF(P,[],isProj);
    [U,S,V]=svd(F); S=diag(S);
    a=sym('a','real');
    w=sym('[a 0 0; 0 a 0; 0 0 1]'); % a is alpha^2
    wp=sym('[a 0 0; 0 a 0; 0 0 1]');
    U=vpa(U); V=vpa(V); S=vpa(S);
    
    eq=cross([U(:,2)'*wp*U(:,2);-U(:,1)'*wp*U(:,2);...
      U(:,1)'*wp*U(:,1)],[S(1)^2*V(:,1)'*w*V(:,1);...
      S(1)*S(2)*V(:,1)'*w*V(:,2);S(2)^2*V(:,2)'*w*V(:,2)]);
    eq=eq(1)^2+eq(2)^2+eq(3)^2;
    coef=sym2poly(eq);
    coefd=sym2poly(diff(eq));
    
    sol=roots(coefd); sol=sol(isreal(sol)); sol=sol(sol>=0);
    [disc,ind]=min(polyval(coef,sol));
    K=sqrt(sol(ind));
  else
    % Reference: HZ2, p163
    P = P/norm(P(3,1:3));
    [ K R ]=rq(P(:,1:3));
    [ K Q ] = makeKPositive(K);
    R = Q*R;
    if nargout>2; t=K\P(:,4); end
  end
else
  % Reference: HZ2, p169
  KR = P(1:2,1:3);
  K = chol( KR*KR' ); K = makeKPositive(K);
  P = K\P(1:2,:);
  R = rotationMatrix(rotationMatrix(P(:,1:3)));
  t = [ P(:,4); 0 ];
  K(3,3)=1;
end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [ K Q ] = makeKPositive(K)
Q = eye(size(K,1));
for i=1:size(K,1); if K(i,i)<0; Q(i,i) = -1; end; end
K = K*Q;
end
