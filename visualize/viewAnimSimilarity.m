function viewAnimSimilarity( anim1, Sim, varargin )
% View the similarities within an anim or between 2 anims
%
% USAGE
%  viewAnimSimilarity( anim1, Sim, varargin )
%
% INPUTS
%  anim     - object of class Animation
%  Sim      - nFrame x nFrame similarity matrix
%  varargin   - list of paramaters in quotes alternating with their values
%   -anim2  - other anim to compare to
%   -nCam   - [-1] Number of cameras to display before and after the
%               current one (-1, don't display the current one)
%   -animGT - ground truth anim
%
% OUTPUTS
%
% EXAMPLE
%
% See also COMPUTEANIMSIMILARITY
%
% Vincent's Structure From Motion Toolbox      Version 3.0
% Copyright (C) 2009 Vincent Rabaud.  [vrabaud-at-cs.ucsd.edu]
% Please email me if you find bugs, or have suggestions or questions!
% Licensed under the GPL [see external/gpl.txt]

global possibleKey showConn showFirst showGT showTitle anim ...
  camMode alignGT hCam cam hPoint hGT hConn;

showGT = false; alignGT = false; camMode = 1;

[ anim2 nCam animGT showFirst ] = getPrmDflt( varargin,{ 'anim2',anim1,...
  'nCam',-1,'animGT',cell(1,2), 'showFirst', false}, 1 );

if ~isa(anim1,'Animation') || ~isa(anim2,'Animation')
  error('anim must be of class Animation');
end

anim = { anim1 anim2 };
cam = cell(1,2); hPoint = cam; hConn = cam; hCam = cam; hGT = cam;

% Determine the boundaries of the data
for i=1:2
  [ cam{i} anim{i} ] = cloudInitializeCam( anim{i}, nCam );
end

% Create the plots
figure(gcf); clf; % bring to focus
h(1) = subplot( 'position', [ 0, 0, 0.3, 0.92 ] );
h(2) = subplot( 'position', [ 0.3, 0, 0.3, 0.92 ] );
h(3) = subplot( 'position', [ 0.6, 0, 0.4, 0.92 ] );
pos=get(gcf,'Position'); pos(1)=pos(1)-1.5*pos(3); pos(3)=3*pos(3);
set(gcf,'Position',pos);

% Show the similarity matrix
imshow( Sim, [] ); hold on;
marker1 = plot( 10, 1, 'r*' ); marker2 = plot( 1, 10, 'b*' );

% Define some initial variables
set( gcf, 'WindowButtonMotionFcn', { @interface } );
set( gcf, 'KeyPressFcn', { @interface } );

possibleKey = { 'f' 'r' 't' '1' '2' '3' 'numpad0' ...
  'leftarrow' 'rightarrow' 'uparrow' 'downarrow' };

hPoint=cell(1,2); hCam=hPoint;

c=[ 1 0.4 0.4; 0.4 0.4 1 ];
for i=1:2
  axes(h(i)); axis on; axis vis3d
  [hPoint{i} hConn{i} hCam{i} hGT{i}]=cloudInitialize( anim{i}, 1, ...
    'nCam',nCam, 'c',c(i,:) );
  for j=1:3; cam{i}(j).gca = h(i); end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  function interface( src, event )
    % Deal with the mouse moving around
    point = get( h(3), 'CurrentPoint' );
    x = round( point( 1, 1:2 ) );
    if any(x<=0) || any(x>[size(Sim,1) size(Sim,2)]); return; end
    
    % Deal with the markers
    set( marker1, 'XData', x(1), 'YData', 1 );
    set( marker2, 'XData', 1, 'YData', x(2) );
    
    for k=1:2
      hCam{k} = cloudUpdate( anim{k}, hPoint{k}, x(k), 'hCam',hCam{k},...
        'nCam', nCam, 'camMode', camMode, ...
        'hGT', hGT{k}, 'hConn',hConn{k},'animGT',animGT{k}, ...
        'showTitle', showTitle, 'showGT', showGT, 'alignGT',...
        alignGT, 'showConn', showConn );
      axis( cam{k}(camMode).gca, cam{k}(camMode).axis );
    end
    
    % Display the image
    cloudInterface( src, event );
  end
end
