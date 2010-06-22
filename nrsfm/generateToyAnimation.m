function anim = generateToyAnimation( type, varargin )
% Generate animation data for SFM
%
% type 0 is a random rigid 3D structure
% type 0.1 is just a simple random linear non-rigid object
% type 1 is a bending cylinder
% type 1.5 is a bending cylinder from blender
% type 2 is a roller coaster
% type 3 is a slinky
% type 3.5 is a real data slinky
% type 4 bending shark from Torresani (NIPS 03)
% type 4.1 sources for the bending shark from Torresani (NIPS 03)
% type 4.2 bending shark from Torresani (NIPS 03) but with random camera
% type 5 is a walking person from Torresani (PAMI 08)
% type 6 is the face from Torresani (PAMI 08)
%
% The camera rotation/translation are defined such that:
%             W = R*S + t
%
% USAGE
%  anim = generateToyAnimation( type, varargin )
%
% INPUTS
%  type       - type of animation to create
%  varargin   - list of paramaters in quotes alternating with their values
%       - 'nPoint' [] number of points to create in the object
%       - 'isProj' [false] flag indicating if the camera is projective
%       - 'nLoop' [1] number of times the animation is repeated
%       - 'dR' [1] set the camera rotation behavior. 0: fixed,
%              1:random (but containing all the scene; if set, dt cannot be
%              set to 2), 2: smooth changes around the scene
%              if a cell containing a matrix, it is supposed to be fixed
%              to that value
%       - 'dt' [1] set the camera translation behavior. Same as for dR
%       - 'nFrame' [10] number of frames in the video sequence
%       - 'K' [eye(3)] calibration matrix
%       - 'randK' [false] set the calibration matrix to random
%       - 'nBasis' [0] dimensionality of the shape basis (for NRSFM)
%           if an array of values, it contains the number of basis shapes
%           of rank 3, 2, 1,   e.g : [ 1 0 1 ]
%       - 'doSMean' [0] if ~=0, Torresani's model, otherwise, Xiao's
%                   (Torresani's : the first shape in the basis always has
%                   a coefficient of 1
%
% OUTPUTS
%  anim           - Animation object
%
% EXAMPLE
%
% See also
%
% Vincent's Structure From Motion Toolbox      Version 3.0
% Copyright (C) 2008-2010 Vincent Rabaud.  [vrabaud-at-cs.ucsd.edu]
% Please email me if you find bugs, or have suggestions or questions!
% Licensed under the GPL [see external/gpl.txt]

dfs = {'nPoint',[], 'isProj',false,'nLoop', 1,...
  'dR',[],'dt',[],'nFrame',10,'K',[],...
  'randK',false,'nBasis',0,'doSMean',0};
[ nPoint isProj nLoop dR dt nFrame K ...
  randK nBasis doSMean ] = getPrmDflt( varargin, dfs, 1 );

% Deal with internal parameters
if randK
  K=eye(3);
  K(1,1) = 500+1000*rand(); K(2,2) = -0.5+rand(); K(1,2)=500+1000*rand();
  if isProj; K(1:2,3) = -0.5 +  rand(); end
end
conn=[];

% Initialize the anim and the camera movements
anim=Animation;
if type>0
  if isempty(dR); dR=2; end; if isempty(dt); dt=0; end
else
  if isempty(dR); dR=1; end; if isempty(dt); dt=0; end
end

switch type
  case 0 % random rigid scene
    if isempty(nPoint) || isempty(nFrame)
      error('Need to specify the number of points and frames');
    end
    S=rand(3,nPoint);
  case 0.1 % random linear non-rigid scene
    if isempty(nPoint) || isempty(nFrame) || isempty(nBasis)
      error('Need to specify the number of points and frames and basis');
    end
    if isProj; error('not implemented yet !'); end
    
    % Create the basis
    SBasis = 2*rand(3,nPoint,sum(nBasis))-1;
    switch length(nBasis)
      case 3,
        SBasis(3,:,nBasis(1)+1:sum(nBasis(1:2)))=0;
        SBasis(2:3,:,sum(nBasis(1:2))+1:end)=0;
      case 1,
      otherwise,
        error('invalid nBasis');
    end
    for i = 1 : size(SBasis,3) % center it
      SBasis(:,:,i) = bsxfun(@minus,SBasis(:,:,i),mean(SBasis(:,:,i),2));
    end
    
    % Create the shape coefficients
    l = zeros( sum(nBasis), nFrame );
    for k = 1 : sum(nBasis)
      %      l(k,:) = smooth( smooth( 2*rand(nFrame,1)-1 ), 25, 'loess' );
      l(k,:)=gaussSmooth( 2*rand(nFrame,1)-1, 6, 'same', 4 );
      l(k,:)=l(k,:)/max(abs(l(k,:)));
    end
    
    if doSMean; l(1,:)=1; end
    anim.l=l; anim.SBasis=SBasis;
    
    % Scale the whole thing
    anim.SBasis=anim.SBasis/max(abs(anim.S(:)))/2;
    S=anim.S;
    
    if dR==1; dR = 2; end
    if dt==1; dt = 0; end
  case 1 % Bending Cylinder
    [ S conn ] = makeCylinder( 8, 5, 200 );
  case 1.5 % Cylinder case
    % read the 3D data
    A = dlmread('./data/cylinder/coord.txt');
    nFrame = nnz(isinf(A(:,2)));
    nPointA = size(A,1)/nFrame - 1;
    
    A(1:nPointA+1:end,:)=[];
    S = reshape(A',3,[],nFrame);
  case 2 % Roller Coaster with a base
    if nFrame==10; nFrame = 150; end
    S = makeCoaster(5,2,nFrame,1);
  case 2.1 % Roller Coaster
    if nFrame==10; nFrame = 150; end
    S = makeCoaster(5,2,nFrame,0);
  case 3 % Slinky
    S = makeSlinky(3,3,500);
  case 3.5 % Real slinky
    load( [ fileparts(mfilename('fullpath')) '/data/slinky/slinky.mat' ] );
    anim.W = normalize( W, nPoint );
    return
  case {4,4.3}, % Shark
    S = makeShark(); dR = 0; dt = 0; isProj = false;
  case 4.1 % Sources of the shark (best result I ever had for 4.1)
    load( [ fileparts(mfilename('fullpath')) '/data/shark/jawSource' ] );
  case 4.2 % Random Camera Shark
    load( [ fileparts(mfilename('fullpath')) '/data/shark/jawSource']);
    
    dR = 2; dt = 0;
    isProj = false;
  case 5 % walking person from Torresani (PAMI 08)
    load( [ fileparts(mfilename('fullpath')) '/data/walkingTorresani08/walking' ] );
    S = permute( reshape( P3_gt, [], 3, size(P3_gt,2) ), [2, 3, 1] );
    dR = 2; dt = 0;
    % add connectivity between points
    conn{1}=[ 26 9 28 26 6 22 21 17 ]; % right leg
    conn{2}=[ 27 25 24 23 ]; % left leg
    conn{3}=[ 1 31 32 33 31 ]; % head
    conn{4}=[ 38 39 43 11 ]; % right arm
    conn{5}=[ 40 41 4 ]; % left arm
    
    for i=1:size(anim.S,2)
      text(anim.S(1,i,1),anim.S(2,i,1),anim.S(3,i,1),int2str(i));
      hold on;
    end
    axis([-1 1 -1 1 -1 1]/2)
  case 6 %face from Torresani (PAMI 08)
  otherwise
    error('type not defined');
end

% Sample
if ~isempty(nPoint) && nPoint<size(S,2); S = S(:,randsample(size(S,2),...
    nPoint),:); end

% Normalize
nFrame=max(nFrame,size(S,3)); nPoint=size(S,2);
if type~=0.1; S = normalize( S, nPoint ); end

% Repeat the anim several times
S=repmat(S,[1 1 nLoop]);
for i=2:2:nLoop
  S(:,:,nFrame*(i-1)+1:nFrame*i) = S(:,:,nFrame:-1:1);
  if ~isempty(anim.l); anim.l(:,nFrame*(i-1)+1:nFrame*i) = anim.l(:,nFrame:-1:1); end
end
nFrame=nLoop*nFrame;

% Create the rotation matrices
if iscell(dR)
  R = repmat(dR{1},[ 1 1 nFrame ]); dR = 0;
else
  R = repmat(eye(3),[ 1 1 nFrame ]);
  switch dR
    case 1, % Purely random
      for i = 1 : nFrame
        vec = rand(1,3); vec = vec/norm(vec)*rand()*2*pi;
        R(:,:,i) = rotationMatrix( vec );
      end
    case 2, % Smooth random
      R = zeros(3,3,nFrame);
      [ R3 uTot ] = circleSpline( nFrame ); % Third rows of the rotation matrices
      R3 = -R3;
      % For each frame, choose one vector ortho to R3 but alined with the previous and next
      % camera positions
      for f = 1 : nFrame
        % find u ortho to R3
        u = uTot(:,f) - (uTot(:,f)'*R3(:,f))*R3(:,f); u = u/norm(u);
        RTmp = [ u cross( R3(:,f), u ) R3(:,f) ]';
        if det(RTmp)<0; RTmp(2,:) = -RTmp(2,:); end
        R(:,:,f) = RTmp;
      end
  end
end

% Create the camera translations
t=zeros(3,nFrame); span = max(abs(S(:)));
if iscell(dt)
  t = repmat(dt{1},[ 1 nFrame ]);
else
  switch dt
    case 1,
      t = [ span*(0.5*rand(2,nFrame)-1); span*(1+rand(1,nFrame)) ];
    case 2,
      t = gaussSmooth( [ span*(0.5*rand(2,nFrame)-1); ...
        span*(1+rand(1,nFrame)) ], 4, 'same', 4 );
  end
end

% make sure all the points are in front of the camera
isGood = false;
while ~isGood
  isGood = true;
  for i = 1 : nFrame
    if size(S,3)==1; STmp=S; else STmp=S(:,:,i); end
    if ~allPointsInFront( STmp, R(:,:,i), t(:,i) ) || ...
        ( ~isempty(S) &&  ...
        ~allPointsInFront( STmp, R(:,:,i), t(:,i) ) )
      % Get the camera a little bit further
      t(3,:) = t(3,:) + span/5; isGood = false;
    end
  end
end

% Define some members
anim.R=R; anim.t=t; anim.K=K;

if isempty(anim.l) || isempty(anim.SBasis); anim.S=S; end

anim.conn=conn; anim.isProj=isProj;

if dR~=0; anim = anim.setFirstPRtToId(); end

[ disc, disc, anim ]=anim.generateW('doFillW',true);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function S = normalize( S, nPoint )
S = S - repmat( mean( S, 2 ), [ 1 nPoint 1 ] );
S = S/max(abs(S(:)))/2;
end
function answer = allPointsInFront( S, R, t )
% compute image points
X = R*S + repmat( t, [ 1 size(S,2) ] );
% Check if they are in front
answer = all( X(3,:) > 0 );
end
