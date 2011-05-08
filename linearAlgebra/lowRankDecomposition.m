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
%  W(rand(m,n)>0.9)=NaN; W=W+rand(m,n)*0.001;
%  [P,Q]=lowRankDecomposition(W,r);
%  norm(P*Q-WFull,'fro')
%
% See also
%
% Vincent's Structure From Motion Toolbox      Version 3.1
% Copyright (C) 2008-2011 Vincent Rabaud.  [vrabaud-at-cs.ucsd.edu]
% Please email me if you find bugs, or have suggestions or questions!
% Licensed under the GPL [see external/gpl.txt]

m=size(W,1); n=size(W,2);
WIsnan=isnan(W);
hasAnyNan=any(WIsnan(:));
nIter=100;
N=zeros(m,nIter*m); NWidth=0;
% figure out the number of top columns we need
nCol=r;
while (nchoosek(nCol,r)) < nIter*m && (nCol<n)
  nCol=nCol+1;
end

% sort the columns by the number of missing elements they have
colCount=sum(WIsnan,1); [disc,colInd]=sort(colCount);
colInd=find(colCount<=colCount(colInd(nCol)));
randInd=1; randArray=randperm(nCol);

for i=1:nIter
  colSubset=colInd(randArray(randInd:randInd+r-1)); randInd=randInd+r;
  if (randInd+r-1)>length(randArray);
    randInd=1; randArray=randperm(nCol);
  end
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
  pMinusR=size(Ai,1)-r; if pMinusR<1; continue; end
  % expand N if needed (for speed issues)
  if NWidth+pMinusR>size(N,2); N(1,NWidth+500)=0; end
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
