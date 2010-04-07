function [hPoint hConn hCam hGT] = cloudInitialize( anim, frameNbr, ...
  varargin )
% Initialize a point cloud animation
%
% USAGE
%  [hPoint hConn hCam hGT]=initializeCloud( varargin )
%
% INPUTS
%  - 'anim' animation object to animate
%  - 'frameNbr' number of the frame to update to
%  varargin   - list of paramaters in quotes alternating with their values
%       - 'nCam' [-1] number of cameras to display in addition to the
%                     current one
%       - 'c' [0.4 0.4 1] color of the features
%       - 'animGT' [] ground truth animation
%       - 'markerScale' [1] scale of the markers
%
% OUTPUTS
%  hPoint     - handle of the animation points
%  hConn      - handle of the connectivity lines of the features
%  hCam       - handle of the cameras
%  hGT        - handle of the ground truth points
%
% EXAMPLE
%
% See also CLOUDUPDATE
%
% Vincent's Structure From Motion Toolbox      Version 2.1
% Copyright (C) 2009 Vincent Rabaud.  [vrabaud-at-cs.ucsd.edu]
% Please email me if you find bugs, or have suggestions or questions!
% Licensed under the Lesser GPL [see external/lgpl.txt]

[ nCam c animGT markerScale ] = getPrmDflt( varargin, ...
  {'nCam',-1,'c',[0.4 0.4 1],'animGT',[],'markerScale',1}, 1 );

hGT=0;
S = anim.S;
cla;

% Initialize the connectivity links
if ~isempty(anim.conn)
  conn=anim.conn;
  hConn = zeros( 1, length(conn) );
  for j=1:length(conn)
    hConn(j)=line( S(1,conn{j},frameNbr), S(2,conn{j},frameNbr), ...
      S(3,conn{j},frameNbr) );
    set(hConn(j),'Visible','off','LineWidth',2*markerScale); hold on;
  end
else
  hConn=0;
end

% Draw the Nodes
if exist('OCTAVE_VERSION','builtin')
  hPoint=plot3( S(1,:,frameNbr), S(2,:,frameNbr), S(3,:,frameNbr), '*', ...
    'Color',c);
else
  hPoint=plot3( S(1,:,frameNbr), S(2,:,frameNbr), S(3,:,frameNbr), '.', ...
    'Color',c );
end
hold on;
hPoint(2)=plot3( S(1,1,frameNbr), S(2,1,frameNbr), S(3,1,frameNbr),'ks');
if ~exist('OCTAVE_VERSION','builtin')
  set(hPoint,'MarkerSize',10*markerScale);
end

if ~isempty(animGT)
  hGT=plot3( animGT.S(1,:,frameNbr),animGT.S(1,:,frameNbr),...
    animGT.S(1,:,frameNbr),'ro');
  hGT(2)=plot3(animGT.S(1,1,frameNbr),animGT.S(1,1,frameNbr),...
    animGT.S(1,1,frameNbr),'rs');
  set(hGT(:),'Visible','off','MarkerSize',5*markerScale);
end

% Initialize the cameras
if nCam>=0
  nFrame = size(anim.t,2);
  nCam = min( nCam, nFrame - 1 );
  hCam=zeros(1+2*nCam,8);
  for j = 1 : 2*nCam + 1
    if exist('OCTAVE_VERSION','builtin')
      for i = 1 : 8
        hCam(j,i) = line( zeros(2,1), zeros(2,1), zeros(2,1) );
      end
    else
      hCam(j,:) = line( zeros(2,8), zeros(2,8), zeros(2,8) );
    end
    set(hCam(j,:), 'Visible', 'off' );
  end
  hCam(:,9) = 0; % Contains the frame they match to
else
  hCam=[];
end

hCam = cloudUpdate( anim, hPoint, frameNbr, 'hCam',hCam,...
  'nCam', nCam, 'hGT', hGT, 'hConn', hConn, 'animGT', animGT );
