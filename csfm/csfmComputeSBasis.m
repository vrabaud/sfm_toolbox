function SBasis = csfmComputeSBasis( l, W, rankS, lamR, addOneToLFlag )
% Compute the best shape matrix in NRSFM given l W and the rank
%
% Just check the Rabaud Belongie CVPR09 paper or Rabaud's thesis
%
% USAGE
%  SBasis = csfmComputeSBasis( l, W, rankS, addOneToLFlag )
%
% INPUTS
%  l             - [ nBasis x nFrame ] shape coefficients
%  W             - [ 2 x nPoint x nFrame ] set of 2D points
%  rankS         - rank of the shape basis
%  addOneToLFlag - if set, 1 is added as the first coefficient of the l
%
% OUTPUTS
%  SBasis        - [ 3 x nPoint x rankS ] shape basis
%
% EXAMPLE
%
% See also
%
% Vincent's Structure From Motion Toolbox      Version 3.0
% Copyright (C) 2008-2011 Vincent Rabaud.  [vrabaud-at-cs.ucsd.edu]
% Please email me if you find bugs, or have suggestions or questions!
% Licensed under the GPL [see external/gpl.txt]

% Predefine some variables
nPoint=size(W,2); nFrame=size(W,3);
if size(W,1)==2; W=reshape(permute(W,[1 3 2 ]),[],nPoint); end
if addOneToLFlag; l=[ ones(1, size(l,2) ); l ]; end

nBasis=size(l,1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
W = bsxfun( @minus, W, mean(W,2) );
[ U S V ] = svd( W, 'econ' );
B=sqrt(S(1:rankS,1:rankS))*V(:,1:rankS)';
WBPlus=W/B;
WtBPlus = permute(reshape(WBPlus',rankS,2,nFrame),[2,1,3]);

lkron = zeros(3, 3*nBasis, nFrame);
for t = 1 : nFrame; lkron(:,:,t) = kron(l(:,t)',eye(3)); end

errBest = Inf;
for globItr = 1 : 20
  G = 10*(rand(3*nBasis, rankS)-0.5);
  
  [ G, R, err ] = optimizeLocal(G,10);
  
  if err<errBest; GBest = G; errBest = err; end
end

% compute a few more iterations with the best G found
[ G, R ] = optimizeLocal(GBest,50);

% Solve for the Q such that Rt*Q is a rotation
QQt = sdpvar(3,3);
T = sdpvar(2*nFrame,2);
F = set('QQt>=0');
for t = 1 : nFrame
  F = F + set('R(:,3*t-2:3*t)*QQt*R(:,3*t-2:3*t)''==eye(2)+T(2*t-1:2*t,:)');
end
obj = norm(T,'fro')^2;%sum(abs(T(:)));
solvesdp( F,obj,sdpsettings('solver','sdpa,csdp,sedumi,*')); % sdpt3

[ U S V ] = svd(double(QQt));
Q = U*real(sqrt(S));

SBasis=permute( reshape( kron(eye(nBasis),inv(Q))*G*B, 3, nBasis, nPoint), [ 1 3 2 ] );


  function [ G, R, err ] = optimizeLocal(G,nItr)
    M = zeros(2*nFrame,3*nBasis);
    
    errPrev = Inf;
    % perform nItr alternate minimization
    for itr = 1 : nItr
      % compute all the Ct
      Ct = multiTimes(lkron,G,1);
      % compute all the Ct*Ct'
      CtCtPrime = multiTimes(Ct,Ct,2.2);
      
      % R AA=BB, let's compute AA
      temp=repmat(lamR*[-1,2,-1],3*nFrame,1);
      AA = spdiags(temp,[-3,0,3],3*nFrame,3*nFrame);
      for i = 1 : 3; AA(i,i) = lamR; AA(end+1-i,end+1-i) = lamR; end
      AA = AA + spBlkDiag(CtCtPrime);
      
      % Let's compute BB
      BB = multiTimes(WtBPlus,Ct,2.2);
      BB = reshape(BB,2,[]);
      
      % Let's finally get R
      R = full(BB/AA);
      
      % and therefore M
      for t = 1 : nFrame
        M(2*t-1:2*t,:) = kron(l(:,t)',R(:,3*t-2:3*t));
      end
      % let's continue and compute G
      G = (M'*M+eye(3*nBasis))\(M'*WBPlus);
      
      err = norm( M*G - WBPlus, 'fro' )^2 + lamR*norm(R(:,1:end-3)-R(:,4:end),'fro')^2 + norm(G,'fro')^2;
      if errPrev-err<0.001*errPrev; break; end
      errPrev = err;
    end
  end

end
