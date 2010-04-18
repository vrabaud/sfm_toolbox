function samp = csfmComputeSample(anim, infimum, nSamp, minDist)
% Compute the triplets to use in CSFM for the GNMDS
%
% Just check the Rabaud Belongie CVPR09 paper or Rabaud's thesis
%
% USAGE
%  samp = csfmComputeSample(anim,infimum,nSamp,minDist)
%
% INPUTS
%  anim          - Animation object
%  infimum       - infimum matrix
%  nSamp         - the number of desired samples
%  minDist       - minimum distance required between frames
%
% OUTPUTS
%  samp          - [ nSamp x 8 ] sample matrix. Each row is like
%                  [ frame1, frame2, frame3, supremum, ...
%                  frame4, frame5, frame6, infimum ]
% EXAMPLE
%
% See also
%
% Vincent's Structure From Motion Toolbox      Version 3.0
% Copyright (C) 2009 Vincent Rabaud.  [vrabaud-at-cs.ucsd.edu]
% Please email me if you find bugs, or have suggestions or questions!
% Licensed under the GPL [see external/gpl.txt]

if nargin<3; minDist = 1; end
nFrame = size(infimum,1);

ticId = ticStatus('n samples',5,1);

samp = zeros(nSamp,8);

tmp = repmat( max(infimum,[],2), [ 1 nFrame ] );
tmpClose = cumsum( exp( 10*( tmp - infimum )./tmp ),2 );

% plot( tmpClose(1,:) ); hold off; imshow(-infimum,[]); hold on ;

count=zeros(nFrame,nFrame);
alpha=1;
for i=1:nSamp
  % Get the 2 triplets of frames
  while 1
    % Get 3 samples that are close to each other
    i1 = mod(i-1,nFrame)+1;
    while 1
      i1 = round( rand()*(nFrame-1) ) + 1;
      s = rand(1,2)*tmpClose(i1,end);
      i2 = find(tmpClose(i1,:)>s(1),1,'first');
      i3 = find(tmpClose(i1,:)>s(2),1,'first');
      
      good=[ i1 i2 i3 ];
      
      if min( [ abs(i1-i2), abs(i1-i3), abs(i2-i3) ] ) < minDist
        continue
      end
      
      % make sure it's not the same pairs called over and over
      if checkCount(count, good, alpha); break; end
    end
    
    % compute the reconstruction of those 3 good pairs
    animBest = computeSMFromW( anim.isProj, anim.W(:,:,[i1 i2 i3]),...
      'method', Inf, 'isCalibrated', true );
    [ err errIndiv ] = animBest.computeError();
    changeCount(count, good, 1);
    
    % deduce the supremum
    errIndiv = (1+0.5)*sqrt(errIndiv(1,:));
    samp(i,5:8) = [ i1 i2 i3 1/3*( (errIndiv(1))^2 + ...
      (errIndiv(2))^2 + ( errIndiv(3))^2 ) ];
    
    for j=1:500
      % Get a bad sample
      i1 = round( rand()*(nFrame-1) ) + 1;
      i2 = round( rand()*(nFrame-1) ) + 1;
      i3 = round( rand()*(nFrame-1) ) + 1;
      
      if min( min( [ abs(i1-i2), abs(i1-i3), abs(i2-i3) ] ) ) < minDist
        continue;
      end
      
      % get the infimum and make sure it is superior to the supremum
      samp(i,1:4) = [ i1 i2 i3 ( infimum(i1,i2) + infimum(i1,i3) + ...
        infimum(i2,i3) )/3 ];
      if samp(i,4) < samp(i,8); continue; end
      
      bad=samp(i,1:3);
      if checkCount(count, bad, alpha); break; end
    end
    if samp(i,4) > samp(i,8); break; end
    changeCount(count, good, -1);
  end
  
  % keep track of the pairs that have been called
  changeCount(count, bad, -1);
  
  %    plot(samp(i,2),samp(i,1),'r*'); plot(samp(i,3),samp(i,1),'r*');
  %    plot(samp(i,6),samp(i,5),'g*'); plot(samp(i,7),samp(i,5),'g*');
  tocStatus( ticId, (i-1)/nSamp );
end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function res = checkCount(count, triplet, alpha)
for ii=1:3
  for jj=ii+1:3
    if count(triplet(ii),triplet(jj))>alpha || ...
        count(triplet(ii),triplet(jj))<-alpha
      res = 0;
      return;
    end
  end
end
res = 1;
end

function changeCount(count, triplet, increment)
for ii=1:3
  for jj=1:3
    count(triplet(ii),triplet(jj)) = count(triplet(ii),triplet(jj)) + ...
      increment;
  end
end
end
