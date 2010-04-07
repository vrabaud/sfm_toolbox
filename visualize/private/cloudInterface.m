function cloudInterface( src, event )
% interface function for several animation display functions
%
% USAGE
%  cloudInterface( src, event )
%
% INPUTS
%
% OUTPUTS
%
% EXAMPLE
%
% See also CLOUDBOUNDARY
%
% Vincent's Structure From Motion Toolbox      Version NEW
% Copyright (C) 2009 Vincent Rabaud.  [vrabaud-at-cs.ucsd.edu]
% Please email me if you find bugs, or have suggestions or questions!
% Licensed under the Lesser GPL [see external/lgpl.txt]

global possibleKey showConn showFirst showGT doReturn showTitle anim ...
  camMode alignGT fps doPause hCam cam hPoint hGT hConn isIndep ...
  showPrettyAxes;

if ~iscell(anim); anim = {anim}; end
if ~iscell(cam); cam = {cam}; end
if ~iscell(hCam); hCam = {hCam}; end
if ~iscell(hPoint); hPoint = {hPoint}; end
if ~iscell(hConn); hConn = {hConn}; end
if ~iscell(hGT); hGT = {hGT}; end

nAnim = numel(anim);
% Deal with a pressed key to change the view or quit the animation
diffView=[0 0];
if ~isfield(event,'Key') || ~any(strcmp(event.Key,possibleKey))
  event.Key='';
end
switch event.Key
  case 'c', % Show the connectivity
    showConn = ~showConn;
  case 'f', % Show the first point
    showFirst = ~showFirst;
  case 'g', % Show the ground truth
    showGT=~showGT;
  case 'i', % make a view independent from the other ones
    isIndep=~isIndep;
  case 'p', % toggle between pretty axes
    showPrettyAxes=~showPrettyAxes;
  case 'q',
    doReturn=1;
  case 'r', % rotate the axes
    for i=1:nAnim
      if ~isempty(anim{i}.S)
        anim{i}.S = anim{i}.S([2 3 1],:,:);
        anim{i}.R = anim{i}.R(:,[2 3 1],:);
        anim{i}.t = anim{i}.t([2 3 1],:);
        for j=1:2
          tmp = [ get(hPoint{i}(j),'XData'); get(hPoint{i}(j),'YData'); ...
            get(hPoint{i}(j),'ZData') ];
          set(hPoint{i}(j),'XData',tmp(2,:));
          set(hPoint{i}(j),'YData',tmp(3,:));
          set(hPoint{i}(j),'ZData',tmp(1,:));
          tmp = [ get(hGT{i}(j),'XData'); get(hGT{i}(j),'YData'); ...
            get(hGT{i}(j),'ZData') ];
          set(hGT{i}(j),'XData',tmp(2,:));
          set(hGT{i}(j),'YData',tmp(3,:));
          set(hGT{i}(j),'ZData',tmp(1,:));
        end
      end
    end
  case 't', % Show the title
    showTitle = ~showTitle;
  case {'1','2','3'} % switch 2D/3D view
    if camMode==str2double(event.Key) && src>=0 && ...
        any(strcmp( 'g', possibleKey ))
      alignGT=~alignGT;
    else
      for i=1:nAnim; cam{i}(camMode).axis = axis; end
      camMode = str2double(event.Key);
      for i=1:nAnim
        switch camMode
          case 1,
            axis( cam{i}(camMode).gca, 'on', 'vis3d' );
            set(gcf,'Color',0.8*[1 1 1]);
            if ~isempty(hCam{i}); set(hCam{i}(:,1:8),'Visible','on'); end
          case 2,
            axis( cam{i}(camMode).gca, 'off', 'equal' );
            set(gcf,'Color',[1 1 1]);
            if ~isempty(hCam{i}); set(hCam{i}(:,1:8),'Visible','off'); end
          case 3,
            axis( cam{i}(camMode).gca, 'on', 'vis3d' );
            set(gcf,'Color',[1 1 1]);
            if ~isempty(hCam{i}); set(hCam{i}(:,1:8),'Visible','on'); end
        end
      end
      
      for i=1:nAnim
        axis( cam{i}(camMode).gca, cam{i}(camMode).axis );
        if camMode~=2
          set( cam{i}(camMode).gca, 'CameraPositionMode', 'auto' );
        else
          set( cam{i}(camMode).gca, 'CameraPosition', [0 -5 0] );
        end
        set( cam{i}(camMode).gca, 'CameraTarget', cam{i}(camMode).target );
      end
    end
  case 'numpad0' % re-center the camera
    for i=1:nAnim; cam{i}(camMode).diff=0; end
  case 'subtract', % slow down the speed
    fps = fps-5;
    if fps<0; fps=1; end
  case 'add', % increase the speed
    fps = fps+5;
  case 'space', % pause
    doPause=~doPause;
  case 'leftarrow',
    diffView = [ 10 0 ];
  case 'rightarrow',
    diffView = [ -10 0];
  case 'uparrow',
    diffView = [ 0 -10];
  case 'downarrow',
    diffView = [ 0 10 ];
end

% Update the look of the data
offon = {'off' 'on'};
for i=1:nAnim
  cam{i}(camMode).diff = cam{i}(camMode).diff + diffView;
  
  % Set some stuff visible or not
  set(hPoint{i}(2),'Visible',offon{1+showFirst});
  set(get(cam{i}(camMode).gca,'Title'),'Visible',offon{1+showTitle});
  
  try
    h = set(hGT{i}(1),'Visible',offon{1+showGT});
    h = set(hGT{i}(2),'Visible',offon{1+showGT*showFirst});
  catch %#ok<CTCH>
  end
  try
    h = set(hConn{i},'Visible',offon{1+showConn});
  catch %#ok<CTCH>
  end
end

if camMode~=2 || isempty(event.Key)
  for i=1:nAnim
    set( cam{i}(camMode).gca, 'View', cam{i}(camMode).view + ...
      cam{i}(camMode).diff);
    
    if showPrettyAxes
      set(cam{i}(camMode).gca,'XTick',[],'YTick',[],'ZTick',[],...
        'LineWidth',3);
    else
      set(cam{i}(camMode).gca,'XTickMode','auto','YTickMode','auto',...
        'ZTickMode','auto','LineWidth',1);
    end
  end
end

% In the case of animation playing, reset the variables to normal
if nAnim==1
  anim = anim{1}; cam = cam{1}; hCam = hCam{1}; hPoint = hPoint{1};
  hConn = hConn{1}; hGT = hGT{1};
end
