function [WHat,coeff]=lowRankDecomposition(W,r)
% Decompose a matrix into the product of two fixed rank matrices (NaN allowed)
%
% W can contain missing entries (containes as NaN)
%
% This is an implementation of:
% Low-Rank Matrix Fitting Based on Subspace Perturbation Analysis
% with Applications to Structure from Motion
% Hongjun Jia, Aleix M. Martinez, PAMI 08
%
% USAGE
%   [P,Q]=lowRankDecomposition(W,r)
%
% INPUTS
%  W     - matrix to decompose. Can contain NaN
%  r     - rank of the decomposition
%
% OUTPUTS
%  WHat     - matrix of rank r such that W=WHat*coeff. Contains no NaN.
%  coeff    - matrix of rank r such that W=WHat*coeff. Contains no NaN.
%
% EXAMPLE
%  m=10; n=5; r=3;
%  W=rand(m,r)*rand(r,n); WFull=W;
%  W(rand(m,n)>0.9)=NaN;
%  [P,Q]=lowRankDecomposition(W,r);
%  norm(P*Q-WFull,'fro')
%
% See also
%
% Vincent's Structure From Motion Toolbox      Version NEW
% Copyright (C) 2008-2010 Vincent Rabaud.  [vrabaud-at-cs.ucsd.edu]
% Please email me if you find bugs, or have suggestions or questions!
% Licensed under the GPL [see external/gpl.txt]

m=size(W,1); n=size(W,2);
WIsnan=isnan(W);
hasAnyNan=any(WIsnan(:));

N=zeros(m,1000); NWidth=0;
for i=1:1000
  colSubset=randperm(n)(1:r);
  Bi=W(:,colSubset);
  % remove rows that have a NaN
  if hasAnyNan
    goodRow=~any(WIsnan(:,colSubset),2);
    Ai=Bi(goodRow,:);
  else
    Ai=Bi;
  end
  % get the null space
  Mi=null(Ai');
  pMinusR=size(Ai,1)-r;
  % expand N if needed (for speed issues)
  if NWidth+pMinusR>size(N,2); N(1,NWidth+100)=0; end
  % add Ni
  if hasAnyNan
    N(goodRow,NWidth+1:NWidth+pMinusR)=Mi(:,1:pMinusR);
  else
    N(:,NWidth+1:NWidth+pMinusR)=Mi(:,1:pMinusR);    
  end
  NWidth=NWidth+pMinusR;
end
N(:,NWidth+1:end)=[];
[U,S,V]=svd(N,'econ');
WHat=U(:,m-r+1:end);

% recover the coefficients
if hasAnyNan
  coeff=zeros(r,n);
  for i=1:n
    col=W(:,i);
    % remove rows that have a NaN
    goodRow=~any(isnan(col),2);
    WHat(goodRow,:)\col(goodRow);
    coeff(:,i)=WHat(goodRow,:)\col(goodRow);
  end
else
  coeff=WHat\W;
end
return

m=10; n=5; r=3;
W=rand(m,r)*rand(r,n); WFull=W;
W(rand(m,n)>0.9)=NaN;
[P,Q]=lowRankDecomposition(W,r);
norm(P*Q-WFull,'fro')
