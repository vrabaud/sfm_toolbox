function S = computeSFromWM( varargin )
% Perform triangulation (3D positions from motion and W)
%
% REFERENCE
%  method 0, isProj==true, is linear triangulation from HZ2, p312
%  method Inf is HZ2, p318, Algorithm 12.1
%
% USAGE
%  S = computeSFromWM( 'W1', W1, 'W2', W2, 'isProj', isProj, 'P1', P1, ...
%       'P2', P2, 'method', method );
%  S = computeSFromWM( 'W', W, 'isProj', isProj, 'P', P, 'method', method )
%
% INPUTS
%  varargin   - list of paramaters in quotes alternating with their values
%       - 'W1',[] [ 2 x nPoint ] 2D projected features in the first image
%       - 'W2',[] [ 2 x nPoint ] 2D projected features in the second image
%       - 'W', [] [ 2 x nPoint x nFrame ] 2D projected features
%       - 'isProj','REQ' flag indicating if the camera is projective
%       - 'P1',[] [3 x 4 ] projection matrix of the first camera
%       - 'P2',[] [3 x 4 ] projection matrix of the second camera
%       - 'P',[] [3 x 4 x nFrame ] projection matrices of the camera
%       - 'K',[eye(3)] [3 x 3 ] calibration matrix
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

dfs = {'W1',[],'W2',[],'W',[],'isProj','REQ','P1',[],'P2',[],'P',[],'K',...
  eye(3),'method',0};
[ W1 W2 W isProj P1 P2 P K method ] = getPrmDflt( varargin, dfs, 1 );

if ~isempty(W1); W = W1; W(:,:,2) = W2; P = P1; P(:,:,2) = P2; end

nFrame = size(W,2); nPoint = size(W,2);

% If uncalibrated or if a third coordinate
if any(any(K~=eye(3))) || size(W,1)==3
  invK = inv(K);
  for i=1:nFrame
    W(1:2,:,i) = normalizePoint( invK*normalizePoint( W(:,:,i), -3 ), 3 );
  end
  W = W( 1:2, :, : );
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
