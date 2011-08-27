function [ Sim QTot ] = computeAnimSimilarity( anim, method, ...
  varargin )
% Compute similarities between frames of a video sequence
%
% This routine computes:
%  - the pair infimum of CSFM if method = -2
%  - the best 2-view reconstruction in MSFM if method = 2
%  - the best 3-view reconstruction in MSFM if method = 3
%
% USAGE
%  [ Sim QTot ] = computeAnimSimilarity( anim, method, Sim )
%
% INPUTS
%  anim     - anim object (see GENERATETOYANIMATION for details)
%  method   - [2] uses pairs of frames, 3, uses triplets
%  varargin   - list of paramaters in quotes alternating with their values
%    for triadic comparison
%       'Sim2' pairwise similarity computed using method = 2
%       'perc' percentage of triplets to consider
%             ( only for speed purposes )
%
% OUTPUTS
%  Sim      - [ nFrame x nFrame ] or [ nFrame x nFrame x nFrame ]similarity
%             matrix based on total reconstruction error
%  QTot     - [ 4 x nFrame x nFrame ] optimal quaternions for pairwise
%
% EXAMPLE
%
% See also GENERATETOYANIMATION, VIEWANIMSIMILARITY
%
% Vincent's Structure From Motion Toolbox      Version NEW
% Copyright (C) 2008-2011 Vincent Rabaud.  [vrabaud-at-cs.ucsd.edu]
% Please email me if you find bugs, or have suggestions or questions!
% Licensed under the GPL [see external/gpl.txt]

if ~isa(anim,'Animation'); error('anim must be of class Animation'); end
nFrame=anim.nFrame;

if nargin<2; method=2; end

% Compute the similarities
switch abs( method )
  case 2
    anim = anim.centerW();

    [ Sim, QTot ] = computeCsfmInfimum(anim.W);
    
    if method==2 % Best 2-view reconstruction for MSFM
      Sim = Sim/2;
    end
    
  case 3
    %%%%%%%% using triplets of frames for MSFM
    [ Sim2 perc ] = getPrmDflt( varargin, { 'Sim2' [] 'perc' 1 } );
    Sim = cell(1,nFrame); for i=1:nFrame; Sim{i}=sparse(nFrame,nFrame); end
    isDone = cell(1,nFrame); for i=1:nFrame; isDone{i}=sparse(nFrame,nFrame); end
    minSpan = 5;
    
    ticId = ticStatus( 'Triadic Computed' );
    
    % compute the probability to choose another frame given a frame in the triplet
    % basically, we'll pick one frame, and the 2 other frames according to that
    % probability. This one simply is proportional to the pairwise reconstruction
    % error (3 frames are more likely to be together if, pairwise, they look alike)
    probS = Sim2;
    
    % remove everything too close to each other
    toep = logical(spdiags(ones(nFrame,2*minSpan+1,-minSpan:minSpan, nFrame,nFrame)));
    probS(toep) = 0;
    
    % compute the cumulative sum and normalize by one
    probS = cumsum(Sim2,2);
    probS = bsxfun(@rdivide,probS,probS(:,end));
    
    % create as many triplets as necessary
    nDone = 0;
    while nDone/nFrame^3*100 < perc
      % choose two other frames according to the distribution
      % in probS and far enough from each other
      sample = randSample(nFrame,3);
      i1 = sample(1); i2 = sample(2); i3 = sample(3);
      if isDone{i1}(i2,i3) || min( pdist( sample' ) )<=minSpan; continue; end
      
      [ disc disc err ] = computeSMFromW( anim.isProj, ...
        anim.W(:,:, [ i1 i2 i3]), 'K', anim.K, 'method',Inf);
      
      perm = [ i1 i2 i3; i1 i3 i2; i3 i1 i2; i3 i2 i1; i2 i1 i3; ...
        i2 i3 i1 ];
      for l=1:6
        Sim{perm(l,1)}(perm(l,2),perm(l,3)) = err(1);
        isDone{perm(l,1)}(perm(l,2),perm(l,3)) = 1;
        nDone = nDone + 1;
      end
      if rand()>0.99
        tocStatus(ticId,min(perc,100*nDone/nFrame^3) / perc);
      end
    end
end
