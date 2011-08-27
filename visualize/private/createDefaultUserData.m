function userData = createDefaultUserData(nAnim)
% create default UserData for the graph
%
% USAGE
%  createDefaultUserData()
%
% INPUTS
%  nAnim  - the number of animations
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

userData.anim = cell(1,nAnim);
userData.cam = cell(1,nAnim);
userData.hCam = cell(1,nAnim);
userData.hPoint = cell(1,nAnim);
userData.hConn = cell(1,nAnim);
userData.hGT = cell(1,nAnim);

userData.showConn = false;
userData.showFirst = false;
userData.showGT = false;
userData.alignGT = false;
userData.isIndep = false;
userData.showPrettyAxes = false;
userData.doReturn = 0;
userData.showTitle = true;
userData.camMode = 1;
userData.fps = 20;
userData.doPause = false;
