function [ cam anim ] = cloudInitializeCam( anim, nCam )
% Initialize the cameras when displaying an animation
%
% USAGE
%  [ cam anim ] = cloudInitializeCam( anim, nCam )
%
% INPUTS
%  anim   - animation object to animate
%  nCam   - [-1] number of cameras to display in addition to the current
%           one. If <0, no camera is displayed
%
% OUTPUTS
%  cam    - cam objects
%  anim   - animation object to animate (with additional data)
%
% EXAMPLE
%
% See also CLOUDBOUNDARY
%
% Vincent's Structure From Motion Toolbox      Version 2.01
% Copyright (C) 2009 Vincent Rabaud.  [vrabaud-at-cs.ucsd.edu]
% Please email me if you find bugs, or have suggestions or questions!
% Licensed under the Lesser GPL [see external/lgpl.txt]

nFrame=anim.nFrame;

if isempty(anim.R), anim.R=repmat(eye(3),[1,1,nFrame]); end
if isempty(anim.t), anim.t=zeros(3,nFrame); end

% Determine the boundaries of the 3D data, camera included
animCam=anim.generateCamFromRt();
cam(1).axis = minmax(reshape(anim.S,3,[]));
if nCam>=0
  cam(1).axis = minmax([cam(1).axis, animCam]);
end

cam=repmat(cam,[1 3]);

% Determine the boundaries of the 2D projected data
ST = anim.generateSAbsolute();

if isempty(anim.W)
  [ disc, disc, anim ]=anim.generateW('doFillW',true);
end

cam(2).axis = [ minmax( reshape(anim.W,2,[]) ); 0 0 ];

try
  cam(3).axis = minmax(reshape(ST,3,[]));
catch %#ok<CTCH>
  cam(3).axis = [];
end

% normalize boundaries
for i=1:3
  maxB=max(cam(i).axis(:,2)-cam(i).axis(:,1))/2;
  % make axes equal
  cam(i).axis=mean(cam(i).axis,2);
  cam(i).axis=[cam(i).axis-maxB cam(i).axis+maxB];
  cam(i).axis=reshape(cam(i).axis',[],1);
end
cam(3).axis(5)=0;
cam(3).axis = cam(3).axis( [ 1 2 5 6 3 4 ] );
cam(2).axis = cam(2).axis( [ 1 2 5 6 3 4 ] );

% define the properties of the different camera modes
cam(1).view = [ -37.5 30 ];
cam(1).diff = [ 0 0 ];
cam(2).view = [ 0 0 ];
cam(2).diff = [ 0 0 ];
cam(3).diff = [ 0 0 ];
cam(3).view = [ 40 20 ];

for i=1:3
  cam(i).target = mean( reshape( cam(i).axis, 2, 3 ), 1 );
end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function mm = minmax( X )
mm = [ min( X, [], 2), max( X, [], 2) ];
end
