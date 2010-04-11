function S = computeSFromWM( isProj, varargin )
% Perform triangulation (3D positions from motion and W)
%
% REFERENCE
%  method 0, isProj==true, is linear triangulation from HZ2, p312
%  method Inf is HZ2, p318, Algorithm 12.1
%
% USAGE
%  S = computeSFromWM( 'W', W, 'isProj', isProj, 'P', P, 'method', method )
%
% INPUTS
%  isProj     - flag indicating if the camera is projective
%  varargin   - list of paramaters in quotes alternating with their values
%       - 'W', [] [ 2 x nPoint x nFrame ] 2D projected features
%       - 'isProj','REQ' flag indicating if the camera is projective
%       - 'P',[] [3 x 4 x nFrame ] projection matrices of the camera
%       - 'K',[eye(3)] [3 x 3 ] or [3 x 3 x nFrame ] calibration matrices
%       - 'method', [0] method for triangulation
%
% OUTPUTS
%  S    - [ 3 x nPoint x nFrame ] 3D structure
%
% EXAMPLE
%
% See also
%
% Vincent's Structure From Motion Toolbox      Version 3.0
% Copyright (C) 2009 Vincent Rabaud.  [vrabaud-at-cs.ucsd.edu]
% Please email me if you find bugs, or have suggestions or questions!
% Licensed under the GPL [see external/gpl.txt]

dfs = {'W',[],'P',[],'K', [],'method',0};
[ W P K method ] = getPrmDflt( varargin, dfs, 1 );

nFrame = size(W,2); nPoint = size(W,2);

% If uncalibrated or if a third coordinate
if ~isempty(K)
  if size(K,3)==1
    invK = repmat( inv(K), [ 1 1 nFrame ] );
  else
    invK = zeros(3,3,nFrame);
    for i=1:nFrame; invK(:,:,i) = inv(K(:,:,i)); end
  end
  W = normalizePoint( multiTimes( invK, normalizePoint( W, -3 ), 2), 3 );
end

switch method
  case 0
    % Linear triangulation
    % Reference: HZ2, p312
    if isProj
      S=zeros(4,nPoint);
      for  i = 1 : nPoint
        % Initial estimate of the 3D points
        A = [ W(1,i,1)*P(3,:,1)-P(1,:,1); W(2,i,1)*P(3,:,1)-P(2,:,1); ...
          W(1,i,2)*P(3,:,2)-P(1,:,2); W(2,i,2)*P(3,:,2)-P(2,:,2)];
        [ disc disc V ] = svd(A,0);
        S( :, i )=V( :, 4 );
      end
      S=normalizePoint(S,4);
    end
  case Inf
    % Reference: HZ2, p318, Algorithm 12.1
end
