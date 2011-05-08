function S = computeSFromWM( isProj, W, P, varargin )
% Perform triangulation (3D positions from motion and W)
%
% REFERENCE
%  method 0, isProj==true, is linear triangulation from HZ2, p312
%  method Inf is HZ2, p318, Algorithm 12.1
%
% USAGE
%  S = computeSFromWM( isProj, W, P, 'method', method )
%
% INPUTS
%  isProj     - flag indicating if the camera is projective
%  W          - [ 2 x nPoint x nFrame ] 2D projected features
%  P          - [3 x 4 x nFrame ] projection matrices of the camera
%  varargin   - list of paramaters in quotes alternating with their values
%       - 'K',[eye(3)] [ 3 x 3 ] or [3 x 3 x nFrame ] calibration matrices
%       - 'method', [0] method for triangulation
%
% OUTPUTS
%  S    - [ 3 x nPoint x nFrame ] 3D structure
%
% EXAMPLE
%
% See also
%
% Vincent's Structure From Motion Toolbox      Version 3.1
% Copyright (C) 2008-2011 Vincent Rabaud.  [vrabaud-at-cs.ucsd.edu]
% Please email me if you find bugs, or have suggestions or questions!
% Licensed under the GPL [see external/gpl.txt]

dfs = {'K',[],'method',0};
[ K method ] = getPrmDflt( varargin, dfs, 1 );

nFrame = size(W,3); nPoint = size(W,2);

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
      A=zeros(2*nFrame,4,nPoint);
      A=[ bsxfun(@minus,bsxfun(@times,permute(W(1,:,:),[3,1,2]),permute(P(3,:,:),[3,2,1])),permute(P(1,:,:),[3,2,1])); bsxfun(@minus,bsxfun(@times,permute(W(2,:,:),[3,1,2]),permute(P(3,:,:),[3,2,1])),permute(P(2,:,:),[3,2,1]))];
      for  i = 1 : nPoint
        % Initial estimate of the 3D points
        [ disc disc V ] = svd(A(:,:,i),0);
        S( :, i )=V( :, 4 );
      end
      S=normalizePoint(S,4);
    end
  case Inf
    % Reference: HZ2, p318, Algorithm 12.1
end
