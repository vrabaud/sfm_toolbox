function hCam = cloudUpdate( anim, hPoint, frameNbr, varargin )
% Update a point cloud animation
%
% USAGE
%  hCam = updateCloud( varargin )
%
% INPUTS
%  - 'anim' animation object to animate
%  - 'hPoint' handle of the points to update
%  - 'frameNbr' number of the frame to update to
%  varargin   - list of paramaters in quotes alternating with their values
%       - 'hCam' [] handle of the cameras
%       - 'nCam' [-1] number of cameras to display in addition to the
%         current one
%       - 'hGT' [0] handle of the ground truth features
%       - 'animGT' [] ground truth animation
%       - 'hConn' [0] handle of the connectivities
%       - 'camMode' [1] camera mode (See PLAYANIM for a description)
%       - 'showGT' [false] show the ground truth points
%       - 'alignGT' [false] align the feature points onto the GT points
%       - 'showConn' [false] show the connectivities
%
% OUTPUTS
%  - hCam    handle of the cameras
%
% EXAMPLE
%
% See also INITIALIZECLOUD
%
% Vincent's Structure From Motion Toolbox      Version 3.0
% Copyright (C) 2008-2010 Vincent Rabaud.  [vrabaud-at-cs.ucsd.edu]
% Please email me if you find bugs, or have suggestions or questions!
% Licensed under the Lesser GPL [see external/lgpl.txt]

dfs = { 'hCam' [] 'nCam' -1 'hGT' 0 'animGT' [] 'hConn' 0 ...
  'camMode' 1 'showGT' false 'alignGT' false 'showConn' false };
[ hCam nCam hGT animGT hConn camMode showGT ...
  alignGT showConn ] = getPrmDflt( varargin, dfs, 1 );

if size(anim.S,3)>1; S = anim.S(:,:,frameNbr); else S=anim.S; end

nPoint = size( anim.S, 2 );

globR = eye(3); globtr = zeros(3,1); globRGT = eye(3); globtrGT=zeros(3,1);

R=anim.R(:,:,frameNbr); t=anim.t(:,frameNbr);

% Align onto the groundtruth
if showGT || alignGT
  RGT=animGT.R(:,:,frameNbr); tGT=animGT.t(:,frameNbr);
  
  ST = R*S + anim.t(:,frameNbr*ones(1,nPoint));
  if size(animGT.S,3)>1; SGT = animGT.S(:,:,frameNbr);
  else SGT=animGT.S;
  end
  
  SGTT = RGT*SGT + repmat( tGT, [1 nPoint] );
  
  if ~alignGT
    [ ST SGTT ] = swap(ST,SGTT); [ t tGT ] = swap(t,tGT);
    [ R RGT ] = swap(R,RGT);
  end
  
  %if alignGT % align on the ground truth
  globR = RGT'*R;
  
  % Chirality ambiguity
  if alignGT
    tempGT = SGTT(3,:) - mean( SGTT(3,:) );
    temp = ST(3,:) - mean( ST(3,:) );
    if norm( temp - tempGT, 'fro') > 1.2*norm( -temp - tempGT,'fro')
      globR = RGT'*[1 0 0; 0 1 0; 0 0 -1]*R;
      t(3) = -t(3);
      ST(3,:) = -ST(3,:);
    end
  end
  temp = t - tGT;
  temp(3) = temp(3) - mean( ST(3,:) ) + mean( SGTT(3,:) );
  
  globtr = RGT'*temp;
  
  if ~alignGT
    [ t tGT ] = swap(t,tGT); [ R RGT ] = swap(R,RGT);
    [ globR globRGT ] = swap(globR,globRGT);
    [ globtr globtrGT ] = swap(globtr,globtrGT);
  end
end

% Align the object to a fixed camera
if camMode==2 || camMode==3
  if alignGT
    globR = RGT*globR; globtr = RGT*globtr + tGT;
    globRGT = RGT; globtrGT = tGT;
  else
    globR = R; globtr = t;
    globRGT = R*globRGT; globtrGT = R*globtrGT + t;
  end
  axisRotation = [ 1 0 0; 0 0 1; 0 1 0];
  globR=axisRotation*globR; globtr=axisRotation*globtr;
  globRGT=axisRotation*globRGT; globtrGT=axisRotation*globtrGT;
end

S = globR*S + repmat( globtr, [ 1 nPoint ] );

% Update the connectivities
if showConn
  conn = anim.conn;
  for j=1:length(hConn)
    set(hConn(j),'XData',S(1,conn{j}),'YData',S(2,conn{j}),...
      'ZData',S(3,conn{j}));
  end
end

% Update the point cloud
set(hPoint(1),'XData',S(1,:),'YData',S(2,:),'ZData',S(3,:));
set(hPoint(2),'XData',S(1,1),'YData',S(2,1),'ZData',S(3,1));

% update the ground truth
if showGT
  SGT = globRGT*SGT + repmat( globtrGT, [ 1 nPoint ] );
  set(hGT(1),'XData',SGT(1,:),'YData',SGT(2,:),'ZData',SGT(3,:));
  set(hGT(2),'XData',SGT(1,1),'YData',SGT(2,1),'ZData',SGT(3,1));
end

% Update the cameras
if nCam>=0
  nFrame = size(anim.t,2);
  nCam = min( nCam, nFrame-1 );
  
  % Redo the frames frameNbr - 1 : frameNbr + 1, no matter what
  % just because of the camera color
  hCam( ismember( hCam(:,9), frameNbr - 1 : frameNbr + 1 ) , 9 ) = 0;

  % Get the frames that are not displayed and that should be displayed
  frameInter = frameNbr - nCam : frameNbr + nCam;
  frameInter(frameInter<1 | frameInter>nFrame) = [];
  toDo = frameInter(~ismember( frameInter, hCam(:,9) ));
  
  % Get the frames whose camera can be deleted
  toErase = find( ~ismember( hCam(:,9), frameInter ) );
  iToErase = 1;
  
  % Process those frames
  for t = toDo
    XX=getCoord(anim.t(:,t),anim.R(:,:,t),0.05+0.05*(t==frameNbr),globR,...
      globtr,camMode);
    nCam = toErase(iToErase); iToErase = iToErase + 1;
    for j=1:8
      set( hCam( nCam, j ),'XData',XX{1}(:,j),'YData',XX{2}(:,j), ...
        'ZData', XX{3}(:,j),'Visible', 'on' );
    end
    if t==frameNbr; set( hCam( nCam, 1 : end-1 ), 'Color', 'r' );
    else set( hCam( nCam, 1 : end-1 ), 'Color', 'b' );
    end
    hCam(nCam,9) = t;
  end
  % Do not display the ones that should not be displayed (if any left)
  hCamToErase=hCam(toErase( iToErase : end ), 1 : end-1);
  set( hCamToErase(:), 'Visible', 'off' );
end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function XX=getCoord(t,R,scale,globR,globtr,camMode)
if (det(globR)<0 && camMode==2) || (det(globR)>0 && camMode==3)
  S=[-1 -1 -2; 1 -1 -2; 1 1 -2; -1 1 -2; 0 0 0]'*scale;
else
  S=[-1 -1 2; 1 -1 2; 1 1 2; -1 1 2; 0 0 0]'*scale;
end
S=(R*globR')'*(S-repmat(t-R*globR'*globtr,[1 5]));

XX=cell(1,3);
for k=1:3
  XX{k}=[S(k,1:2); S(k,2:3); S(k,3:4); S(k,[4 1]); S(k,[1 5]); ...
    S(k,[2 5]); S(k,[3 5]); S(k,[4 5])]';
end
end
function [x y] = swap(x,y)
temp=x; x=y; y=temp;
end
