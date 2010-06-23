function varargout = montageView( S, varargin )
% Used to display collections of 2D/3D points
%
% The following keys can be used when playing an animtion:
%  'c' toggle the connectivity display
%  'f' toggle the first point display
%  't' toggle the title display
%  'numpad0' set the view back to the original
%  {'leftarrow', 'rightarrow', 'uparrow', 'downarrow'} rotate the 3D view
%  'i' toggle the independence of the current view changes
%
% USAGE
%   montageView( S, varargin )
%   [ h m n ] = montageView( S, varargin )
%
% INPUTS
%  S          - [ 2 or 3 x nPoint x nView ]
%  varargin   - list of paramaters in quotes alternating with their values
%       - 'm' [] number of image columns
%       - 'n' [] number of image rowss
%       - 'label' [] cell array of label (strings)
%       - 'conn' [] connecivity array (see GENERATETOYANIMATION)
%
% OUTPUTS
%  h           - plot handle
%  m           - number of image columns
%  n           - number of image rows
%
% EXAMPLE
%
% See also
%
% Vincent's Structure From Motion Toolbox      Version 2.1
% Copyright (C) 2009 Vincent Rabaud.  [vrabaud-at-cs.ucsd.edu]
% Please email me if you find bugs, or have suggestions or questions!
% Licensed under the GPL [see external/gpl.txt]

[ m n label conn ] = getPrmDflt( varargin, { 'm', [], 'n', [], ...
  'label', [], 'conn', [] }, 1 );

nFrame=size(S,3); camMode = 1; showPrettyAxes=false;

% get layout of plots (m=#frames/row, n=#frames/col)
if( isempty(m) || isempty(n))
  if( isempty(m) && isempty(n))
    n = min( ceil(sqrt(nFrame)), nFrame ); m = ceil( nFrame/n );
  elseif( isempty(m) )
    n = min( n, nFrame ); m = ceil(nFrame/n);
  else
    m = min( m, nFrame ); n = ceil(nFrame/m);
  end
end

% Plot the 2D/3D data
k=1; h=zeros(1,nFrame);
hPoint = cell( 1, nFrame ); hConn = hPoint; animTot = hPoint; cam = hPoint;
for j=1:m
  for i=1:n
    h(k)=subplot('position',[(i-1)/n,1-1/m-(j-1)/m,1/n,1/m]);
    if k<=nFrame
      anim=Animation; anim.S=S(:,:,k);
      [ cam{k} animTot{k} ] = cloudInitializeCam( anim, -1 );
      cam{k}(1).gca = gca;
      animTot{k}.conn = conn;
      [ hPoint{k} hConn{k} ] = cloudInitialize( animTot{k}, 1 );
      
      axis( gca, cam{k}(1).axis );
      if ~isempty(label); title(label{k}); end
      k=k+1;
    else % cross out unused frames
      set(gca,'Visible','off'); line( [0 1], [0 1] ); line( [0 1],[1 0] );
    end
  end
end

% Make the whole plot interactive
if size(S,1)==3
  set( gcf, 'WindowButtonMotionFcn', { @interface } );
  set( gcf, 'KeyPressFcn', { @interface } );
end
isIndep=false; showFirst = false; showTitle = false;

possibleKey = { 'f' 'p' 'r' 't' 'numpad0' 'leftarrow' 'rightarrow' ...
  'uparrow' 'downarrow' };
if isempty(hConn{1}); showConn = false; else possibleKey(end+1) = {'c'};end

% optional output
if( nargout>0 ); varargout={ h m n }; end

set(gcf,'UserData',{possibleKey,showConn,showFirst,showTitle,...
  showPrettyAxes,anim,camMode,camm,hPoint,hConn,isIndep});

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  function interface( src, event )
    tmp=set(gcf,'UserData');
    [possibleKey,showConn,showFirst,showTitle,...
      showPrettyAxes,anim,camMode,camm,hPoint,hConn,isIndep]=tmp(:);

    if ~isIndep; inter=1:nFrame; else inter=find(h==gca); end
    
    anim = animTot( inter );
    
    if isempty(event)
      % Set all the views to be like the current one when using the mouse
      cl = find(h==gca);
      if ~isempty(cl) && cl<=nFrame;
        for l=1:nFrame; cam{l}(1).view = get( h(cl), 'View' ); end
        cloudInterface( src, event );
      end
    end
    cloudInterface( src, event );
  end
end
