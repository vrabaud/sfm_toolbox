function [ l K ] = csfmGnmds(samp,infimum,nTriplet,nPair,lam1 )
% Compute the pairs and triplets to use in SDP in CSFM
%
% Just check the Rabaud Belongie CVPR09 paper or Rabaud's thesis
%
% USAGE
%  samp = csfmComputeSample(anim,infimum,nSamp,minDist)
%
% INPUTS
%  anim          - Animation object
%  infimum       - infimum matrix
%  nTriplet      - the number of triplets to use in GNMDS
%  nPair         - the number of pairs to use in GNMDS
%  lam1          - regularization parameter for the smoothness
%
% OUTPUTS
%  l             - [ nBasis x nFrame ] shape coeficients
%  K             - [ nFrame x nFrame ] full rank solution of GNMDS
%
% EXAMPLE
%
% See also
%
% Vincent's Structure From Motion Toolbox      Version 3.0
% Copyright (C) 2009 Vincent Rabaud.  [vrabaud-at-cs.ucsd.edu]
% Please email me if you find bugs, or have suggestions or questions!
% Licensed under the GPL [see external/gpl.txt]

t=find(samp(:,1)==0,1,'first');
if ~isempty(t); samp=samp(1:t-1,:); end

nFrame = size(infimum,1);

% Create Samples
nTriplet = min( [ 3000 size(samp,1) nTriplet ] );
nPair=ceil(sqrt(nPair))^2;
nSampPer3=5;

t=vect(samp(:,[1:3,5:7]),'v');

%%% Start defining the problem
K=sdpvar(nFrame, nFrame);
F=set('K>=0','sd');
obj=0;

%%% Deal with different constraints
% Deal with triplets
if nTriplet>0
  AGood=sparse( nTriplet, nFrame^2 ); ABad=sparse( nTriplet, nFrame^2 );
  for i=1:nTriplet
    indGood=[ index( samp(i,5), samp(i,6) ) index( samp(i,6), samp(i,7) ) index( samp(i,7), samp(i,5) ); ...
      index( samp(i,5), samp(i,5) ) index( samp(i,6), samp(i,6) ) index( samp(i,7), samp(i,7) ) ];
    indBad=[ index( samp(i,1), samp(i,2) ) index( samp(i,2), samp(i,3) ) index( samp(i,1), samp(i,3) ); ...
      index( samp(i,1), samp(i,1) ) index( samp(i,2), samp(i,2) ) index( samp(i,3), samp(i,3) ) ];
    
    tmp=[-1 1];
    for j=1:2
      AGood(i,indGood(j,:) )=AGood(i,indGood(j,:))+tmp(j);
      ABad(i,indBad(j,:) )=ABad(i,indBad(j,:))+tmp(j);
    end
  end
  x1=sdpvar(nTriplet, 1);
  
  F=[ F set('ABad*K(:) + x1 >= AGood*K(:)','triplet') set('x1>=0', 'triplet') ];
  obj=1/nTriplet*sum(x1);
end

% Deal with the minima of pairs
if nPair>0
  A2=sparse( nPair, nFrame^2 );
  b2=zeros( nPair, 1 );
  
  % make sure the maxima are dealt with
  pair = nonMaxSupr( infimum, 1, [], nPair );
  
  % Deal with a sparse grid over infimum
  step=max(floor(sqrt((nPair-size(pair,1))/2)),2);
  [ i, j ] = meshgrid( round(1 : (nFrame-1)/(step-1) : nFrame), ...
    round(1 : (nFrame-1)/(step-1) : nFrame) );
  tmp = [ i(:), j(:) ];
  tmp(tmp(:,1)==tmp(:,2),:) = [];
  pair = [ pair; tmp ];
  
  % finish with random samples
  indStart = size(pair,1)+1;
  if indStart<nPair; pair(nPair,:) = 0; end
  for ind = indStart : nPair
    pair(ind,:) = randsample(nFrame,2);
  end
  
  for i = 1 : nPair
    j = pair(i,:);
    A2(i, [ index( j(1), j(1) ) index( j(2), j(2) )  ...
      index( j(1), j(2) ) ] ) = [ 1 1 -2 ];
    b2(i)=infimum( j(1), j(2) );
  end
  
  F=[ F set('A2*K(:) >= b2 ','pair') ];
end

% centering constraint
F=[ F set('ones(1,nFrame^2)*(K(:)) == 0','center') ];

% Low rank/trace constraint, does not help at all
%  obj=obj+trace(K);

%%% Define the rest of the objective function
% Temporal smoothness regularization

cSmooth1=sparse(1,nFrame^2);
for i = 2 : nFrame
  tmp=[ index( i, i ) index( i, i-1 ) index( i-1, i-1 ) ];
  cSmooth1(tmp)=cSmooth1(tmp)+[ 1 -2 1 ]*lam1;
end
obj=obj+cSmooth1*K(:);

% Temporal smoothness regularization, second order (not necessary)
if 0
  cSmooth2=sparse(1,nFrame^2);
  for i = 2 : nFrame-1
    tmp = [ index( i+1, i+1 ), index( i+1, i ), index( i+1, i-1 ), ...
      index( i, i+1 ), index( i, i ), index( i, i-1 ), ....
      index( i-1, i+1 ), index( i-1, i ), index( i-1, i-1 ) ];
    cSmooth2(tmp)=cSmooth2(tmp) + [ 1 -2 1 -2 4 -2 1 -2 1 ]*lam1;
  end
  obj=obj+cSmooth2*K(:);
end

%%% Solve the whole problem
solvesdp( F,obj,sdpsettings('solver','sdpa,csdp,sedumi,*',...
  'dualize',1,'debug',1));

K = double(K);

[ U S V ] =svd(full(K), 'econ');
l = sqrt(S)*V';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  function indAB = index( a, b )
    indAB = sub2ind( [nFrame nFrame], a, b );
  end
end
