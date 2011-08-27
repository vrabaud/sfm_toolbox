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
% Copyright (C) 2008-2010 Vincent Rabaud.  [vrabaud-at-cs.ucsd.edu]
% Please email me if you find bugs, or have suggestions or questions!
% Licensed under the Lesser GPL [see external/lgpl.txt]

userData=get(gcf,'UserData');

anim = userData.anim; cam = userData.cam; hCam = userData.hCam;
hPoint = userData.hPoint; hConn = userData.hConn; hGT = userData.hGT;

if userData.isIndep
  animRange = 1;
else
  animRange = 1:numel(anim);
end

if ~iscell(hGT); hGT = {hGT}; end

% Deal with a pressed key to change the view or quit the animation
diffView=[0 0];
if ~isfield(event,'Key') || ~any(strcmp(event.Key,userData.possibleKey))
  event.Key='';
end
switch event.Key
  case 'c', % Show the connectivity
    userData.showConn = ~userData.showConn;
  case 'f', % Show the first point
    userData.showFirst = ~userData.showFirst;
  case 'g', % Show the ground truth
    userData.showGT = ~userData.showGT;
  case 'i', % make a view independent from the other ones
    userData.isIndep=~userData.isIndep;
  case 'p', % toggle between pretty axes
    userData.showPrettyAxes=~userData.showPrettyAxes;
  case 'q',
    userData.doReturn=1;
  case 'r', % rotate the axes
    for i=animRange
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
    userData.showTitle = ~userData.showTitle;
  case {'1','2','3'} % switch 2D/3D view
    if userData.camMode==str2double(event.Key) && src>=0 && ...
        any(strcmp( 'g', userData.possibleKey ))
      userData.alignGT = ~userData.alignGT;
    else
      for i=animRange; cam{i}(userData.camMode).axis = axis; end
      userData.camMode = str2double(event.Key);
      for i=animRange
        switch userData.camMode
          case 1,
            axis( cam{i}(userData.camMode).gca, 'on', 'vis3d' );
            set(gcf,'Color',0.8*[1 1 1]);
            if ~isempty(hCam{i})
              set(reshape(hCam{i}(:,1:8),[],1),'Visible','on');
            end
          case 2,
            axis( cam{i}(userData.camMode).gca, 'off', 'equal' );
            set(gcf,'Color',[1 1 1]);
            if ~isempty(hCam{i})
              set(reshape(hCam{i}(:,1:8),[],1),'Visible','off');
            end
          case 3,
            axis( cam{i}(userData.camMode).gca, 'on', 'vis3d' );
            set(gcf,'Color',[1 1 1]);
            if ~isempty(hCam{i})
              set(reshape(hCam{i}(:,1:8),[],1),'Visible','on');
            end
        end
      end
      
      for i=animRange
        axis( cam{i}(userData.camMode).gca, cam{i}(userData.camMode).axis);
        if userData.camMode~=2
          set( cam{i}(userData.camMode).gca, 'CameraPositionMode', 'auto');
        else
          set( cam{i}(userData.camMode).gca, 'CameraPosition', [0 -5 0] );
        end
        set( cam{i}(userData.camMode).gca, 'CameraTarget', ...
          cam{i}(userData.camMode).target );
      end
    end
  case 'numpad0' % re-center the camera
    for i=animRange; cam{i}(userData.camMode).diff=0; end
  case 'subtract', % slow down the speed
    userData.fps = userData.fps-5;
    if userData.fps<0; userData.fps=1; end
  case 'add', % increase the speed
    userData.fps = userData.fps+5;
  case 'space', % pause
    userData.doPause = ~userData.doPause;
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
for i=animRange
  cam{i}(userData.camMode).diff = cam{i}(userData.camMode).diff + diffView;
  
  % Set some stuff visible or not
  set(hPoint{i}(2),'Visible',offon{1+userData.showFirst});
  set(get(cam{i}(userData.camMode).gca,'Title'),'Visible',...
    offon{1+userData.showTitle});
  
  if (~isempty(hGT)) && (~isempty(hGT{i}))
    if hGT{i}(1)>=0
      set(hGT{i}(1),'Visible',offon{1+userData.showGT});
    end
    if (length(hGT{i}) >=2 ) && (hGT{i}(2)>=0)
      set(hGT{i}(2),'Visible',offon{1+userData.showGT*userData.showFirst});
    end
  end
  if (~isempty(hConn))
    set(hConn{i},'Visible',offon{1+userData.showConn});
  end
end

if userData.camMode~=2 || isempty(event.Key)
  for i=animRange
    set( cam{i}(userData.camMode).gca, 'View', ...
      cam{i}(userData.camMode).view + cam{i}(userData.camMode).diff);
    
    if userData.showPrettyAxes
      set(cam{i}(userData.camMode).gca,'XTick',[],'YTick',[],'ZTick',[],...
        'LineWidth',3);
    else
      set(cam{i}(userData.camMode).gca,'XTickMode','auto','YTickMode',...
        'auto','ZTickMode','auto','LineWidth',1);
    end
  end
end

% Save the data back
userData.anim = anim; userData.cam = cam; userData.hPoint = hPoint;
userData.hGT = hGT;

set(gcf,'UserData',userData);
