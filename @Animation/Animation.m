function anim=Animation(varargin)
% Animation class (with 3D data, projection, camera info ...)
%
% An Animation object has the following properties
%   W        - [ 2 x nPoint x nFrame ] matrix containing the projected
%              animation
%   K        - calibration matrix parameters. Can be [ 3 x nFrame ] or
%              [ 5 x nFrame ] or empty. If all the parameters are the same,
%              use [ 3 x 1 ] and [ 5 x 1 ]. The matrix lists the
%              calibration parameters in the following order
%              Affine
%               K1  K2  0
%               0   K3  0
%               0   0   1
%              Projective
%               K1  K2  K4
%               0   K3  K5
%               0   0   1
%   S        - [ 3 x nPoint x nFrame ] matrix containing the 3D
%              animation
%              or [ 2 x nPoint x nFrame ] if working with homographies
%   P        - [ 3 x 4 x nFrame ] projection matrices
%              (can be [ 3 x 3 x nFrame ] in the case of homographies)
%              If it is empty, only K,R,t are used
%   R        - [ 3 x 3 x nFrame ] camera rotation. Can be empty (only P is
%              used)
%   t        - [ 3 x nFrame ] camera translation. Can be empty (only P is
%              used)
%   cam      - [ 3 x nFrame ] position of the camera
%   mask     - [ nPoint x nFrame ] mask of appearance: 1 the point appears
%              in the frame, 0 it does not
%   conn     - cell of arrays indicating connectivities: each point array
%              is a list of points forming a broken line ([1 2 3] means
%              there will be a line from 1 to 2, one from 2 to 3.
%              only fordisplay purposes
%   l        - [ dimSSpace x nFrame ] linear coefficients of the shape
%              basis in NRSFM
%   SBasis   - [ 3 x nPoint x dimSSpace ] Xiao's shape basis or
%              [ 3 x nPoint x (dimSSpace+1) ] Torresani's shape basis
%              (first shape has 1 as first coefficient, and the first
%              coefficient does not need to be specified in l)
%   misc     - whatever you want :)
%   type     - 0: rigid
%              1: rigid with homographies (S is [ 2 x nPoint x nFrame ])
%              2: non-rigid
%   nBasis   - number of shape bases
%   nPoint   - number of feature points
%   nFrame   - number of frames
%   isProj   - flag indicating if the camera is projective
%
%
%
% IMPORTANT
%   for coding comfort, the automatic processes happens:
%     - when l or SBasis is specified, S cannot be changed by an outsider
%     - when any element in l or SBasis is modified, the corresponding
%       elements in S are modified
%
% An Animation object has the following methods
%   anim=generateSFromLSBasis( anim );
%   anim=generateCamFromRt( anim );
%   anim=generateWFromSRt( anim );
%   anim=generateSPCA( anim, nPCA );
%   anim=setFirstRToId( anim );
%   [ anim meanW ]=centerW( anim );
%   [ anim SimMin ]=sampleFrame( anim, nSamp, SimMin );
%
% USAGE
%  anim = Animation()
%
% INPUTS
%   nFrame  - the number of frames to consider
%   nPoint  - the number of points to consider
%   nBasis  - the number of shape bases if any
%   nDim    - the dimensionality of the points (3 for 3D, 2 for 2D)
%
% OUTPUTS
%  anim     - an Animation object
%
% EXAMPLE
%
% See also GENERATETOYANIMATION
%
% Vincent's Structure From Motion Toolbox      Version NEW
% Copyright (C) 2009 Vincent Rabaud.  [vrabaud-at-cs.ucsd.edu]
% Please email me if you find bugs, or have suggestions or questions!
% Licensed under the Lesser GPL [see external/lgpl.txt]

[ nFrame nPoint nBasis nDim ] = getPrmDflt( varargin,{ 'nFrame', 0, ...
  'nPoint', 0, 'nBasis', 0, 'nDim', 0 }, 1 );

% measurements
anim.W=zeros(2,nPoint,nFrame);
anim.mask= false(nPoint,nFrame);

% 3D / Object
anim.S=[];
anim.conn={}; % Cell of arrays of connectivities

% NRSFM
if nDim==3 && nBasis>0
  anim.l=zeros(nBasis,nFrame); anim.SBasis=zeros(nDim,nPoint,nBasis);
else anim.l=[];anim.SBasis=[];
end

% Camera
anim.P=[];
anim.K=[];
if nDim==3; anim.R=zeros(3,3,nFrame); anim.t=zeros(3,nFrame);
else anim.R=[]; anim.t=[]; end
anim.cam=[]; % Position of the optical center

% Misc
anim.misc=[];

% Info
anim.isProj=false;
anim.type=2;

% Quantities
anim.nBasis=nBasis;
anim.nPoint=nPoint;
anim.nFrame=nFrame;

anim = class( anim, 'Animation');
