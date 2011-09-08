function M = playAnim( animIn, varargin )
% Display a point cloud animation
%
% The following keys can be used when playing an animtion:
%  'c' toggle the connectivity display
%  'f' toggle the first point display
%  'g' toggle the ground truth display
%  'p' set the axes to be pretty
%  'q' quit the playing of the animation
%  'r' rotate the axis order
%  't' toggle the title display
%  {'1','2','3'} switch the camera view (1: fixed axes, 2: camera
%                projection, 3:fixed camera}
%  'numpad0' set the view back to the original
%  '-' slow down the speed
%  '+' increase the speed
%  'space' pause the animation
%  {'leftarrow', 'rightarrow', 'uparrow', 'downarrow'} rotate the 3D view
%
% USAGE
%   playAnim( animIn, varargin )
%   M = playAnim( anim, varargin )
%
% INPUTS
%  anim     - object of class Animation (help Animation for details)
%  varargin  - list of paramaters in quotes alternating with their values
%       - 'nCam' [0] number of cameras to show before and after the current
%         one (-1 => not even the current one is displayed)
%       - 'fps' [20] maximum number of frames to display per second
%         use fps==0 to introduce no pause and have the movie play as
%         fast as possible
%       - 'nLoop' [1] number of time to nLoop video (may be inf),
%         if neg plays video forward then backward then forward etc.
%       - 'animGT' [] ground truth animation
%       - 'frame' [] only display the specified frames
%       - 'camMode' [1] choose the original camera mode as camMode
%       - 'showFirst' [false] display the first point differently (with a
%         square)
%       - 'showTitle'  [true] display the title (frame number)
%       - 'showGT' [false] display the ground truth from the very beginning
%       - 'alignGT' [false] align onto the ground truth animation from the
%         very beginning
%       - 'showConn' [true] show the connectivities in anim
%       - 'showPrettyAxes' [false] show pretty axes
%
% OUTPUTS
%   M      - MATLAB movie of the animation
%
% EXAMPLE
%
% See also
%
% Vincent's Structure From Motion Toolbox      Version 3.1.1
% Copyright (C) 2008-2011 Vincent Rabaud.  [vrabaud-at-cs.ucsd.edu]
% Please email me if you find bugs, or have suggestions or questions!
% Licensed under the GPL [see external/gpl.txt]

if ~strcmp( class(animIn), 'Animation' );
  error(' Need to input an Animation object');
end

anim=animIn;
nPoint=anim.nPoint; nFrame=anim.nFrame;
userData = createDefaultUserData(1);

dfs = {'nCam',0,'fps',20, 'nLoop',1, 'animGT',[],'frame',[],'camMode',1,...
  'showFirst',false,'showTitle',true,'showGT',false,'alignGT',false,...
  'showConn',true, 'showPrettyAxes', false };
[ nCam userData.fps nLoop animGT frame userData.camMode ...
  userData.showFirst userData.showTitle userData.showGT userData.alignGT...
  userData.showConn userData.showPrettyAxes ] = getPrmDflt(varargin,dfs,1);

if ~isempty(frame); userData.showGT = true; end
if isempty(anim.conn); userData.showConn=false; end
if isempty(anim.S);
  userData.camMode = 2;
  anim.S = anim.W;
  anim.S(3,:,:) = anim.W(2,:,:);
  anim.S(3,:,:)=0;
  nFrame=size(anim.S,3);
end

canGT = ~isempty(animGT);

if isempty(frame); frame = 1:nFrame; end

if canGT
  if isempty(animGT.S)
    for i=1:nFrame
      animGT.S(:,:,i) = anim.R(:,:,i)'*( [ animGT.W(:,:,i);
        ones(1,nPoint) ] - repmat( anim.t(:,i), [ 1 nPoint ] ) );
    end
    animGT.R = anim.R; animGT.t = anim.t;
  end
end

% Determine the boundaries of the data
[ cam anim ] = cloudInitializeCam( anim, nCam );

% Define some initial variables
h=gcf; figure(h); clf;
set( gcf, 'KeyPressFcn', {  } );
userData.doReturn=0; userData.doPause=0;

userData.possibleKey = { 'f' 'p' 'q' 'r' 't' '1' '2' '3' 'numpad0' ...
  'subtract' 'add' 'space' 'leftarrow' 'rightarrow' 'uparrow' 'downarrow'};

[ userData.hPoint{1} userData.hConn{1} userData.hCam{1} ...
  userData.hGT{1} ]=cloudInitialize( anim, frame(1), 'nCam',nCam, ...
  'c', [0.4,0.4,1], 'animGT', animGT,'markerScale',1 );

% set the main handle
set(gcf,'KeyPressFcn',{@cloudInterface});
if userData.hConn{1}==0; userData.showConn = false;
else userData.possibleKey(end+1) = {'c'};
end
if userData.hGT{1}==0; userData.showGT = false;
else userData.possibleKey(end+1) = {'g'};
end

axis(cam(userData.camMode).axis);
for i=1:3; cam(i).gca = gca; end
% set the data embedded in the figure
userData.anim{1}=anim; userData.cam{1}=cam;
set(gcf,'UserData',userData);
cloudInterface(-1,struct('Key',int2str(userData.camMode)));
cloudInterface(-1,struct('Key','' ));

% play the animation several times
if nargout>0; f=1; M = repmat( getframe, 1, abs(nLoop)*length(frame) ); end
for iLoop = 1 : abs(nLoop)
  % Play the animation once
  for frameNbr=frame
    tic
    while 1
      % embed the right data in the graph
      userData = get(gcf,'UserData');
      try geth=get(h); catch; return; end %#ok<CTCH,NASGU>
      if userData.doReturn; return; end
      
      if userData.showTitle
        set(get(gca,'Title'),'String',...
          sprintf('frame %d of %d',frameNbr,nFrame) );
      end
      
      % Display the image
      userData.hCam = cloudUpdate( userData.anim{1}, userData.hPoint{1}, ...
        frameNbr, 'hCam',userData.hCam{1},'nCam', nCam, 'camMode', ...
        userData.camMode, 'hGT', userData.hGT{1}, 'hConn', ...
        userData.hConn{1}, 'animGT', animGT, 'showGT', userData.showGT, ...
        'alignGT',userData.alignGT, 'showConn', userData.showConn );
      
      if userData.doPause; pause(0.01); continue; end
      
      if  toc>1/userData.fps; break;
      else
        if 1/userData.fps-toc>0.01; pause(0.01);
        else pause(1/userData.fps-toc); break;
        end
      end
    end
    
    if nargout>0; M(f) = getframe(gcf); f=f+1; else drawnow; end
  end
  frame = frame( end:-1:1 );
end
