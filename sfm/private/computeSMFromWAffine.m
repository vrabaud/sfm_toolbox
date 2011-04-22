function [P,S] = computeSMFromWAffine( W, method, onlyErrorFlag )
% Compute SfM from measurements for affine cameras
%
% This function should not be called directly: it does not contain
% pre-processing steps (removal of calibration matrix if known) or
% post-procssing steps (affine/metric upgrades) and bundle adjustment
% computeSMFromW should always be used instead.
%
% Gold Standard for Affine camera matrix
% Reference: HZ2, p351, Algorithm 14.1
% nFrame==2 && onlyErrorFlag
% fast
%
% Tomasi Kanade without the metric constraint
% Affine camera matrix, MLE estimation (Tomasi Kanade)
% Reference: HZ2, p437, Algorithm 18.1
% nFrame>2
% fast
%
% If there are any missing entries:
% Low-Rank Matrix Fitting Based on Subspace Perturbation Analysis
% with Applications to Structure from Motion
% Hongjun Jia, Aleix M. Martinez, PAMI 08
%
% Returns the structure and camera parameters in an Animation object
%
% USAGE
%   anim = computeSMFromWProjective( W, method )
%
% INPUTS
%  W             - [ 2 x nPoint x nFrame ] 2D projected features
%  method        - [inf] method for performing SFM (see above for details)
%  onlyErrorFlag - flag indicating if we only want the error (only used
%                 for 2 frames)
%
% OUTPUTS
%  P          - [ 3 x 4 x nFrame ] projection matrices
%  S          - [ 3 x nPoint ] structure matrix
%
% EXAMPLE
%
% See also
%
% Vincent's Structure From Motion Toolbox      Version NEW
% Copyright (C) 2008-2010 Vincent Rabaud.  [vrabaud-at-cs.ucsd.edu]
% Please email me if you find bugs, or have suggestions or questions!
% Licensed under the GPL [see external/gpl.txt]

nFrame = size(W,3); nPoint = size(W,2);

P=zeros(3,4,nFrame); P(3,3:4)=[0 1];

WIsnan=isnan(W);
hasAnyNan=any(WIsnan(:));

% Gold Standard for Affine camera matrix
% Reference: HZ2, p351, Algorithm 14.1
% Just the error is computed
if nFrame==2 && onlyErrorFlag
  % Normalize input data
  A=[ W(:,:,2); W(:,:,1) ]';
  A=bsxfun(@minus,A,mean(A,1));
  
  N = solveLeastSqAx( A, [], 1 );
  errFrame = N'*A';
  P = norm(errFrame)^2;
else
  WStack = reshape( permute( W, [ 1 3 2 ] ), [], nPoint );
  
  if hasAnyNan
    % Low-Rank Matrix Fitting Based on Subspace Perturbation Analysis
    % with Applications to Structure from Motion
    % Hongjun Jia, Aleix M. Martinez, PAMI 08
    [PStack,S]=lowRankDecomposition(WStack,4);
    % Reconstruct the full measurements
    WStack=PStack*S;
  end
  
  t = mean( WStack, 2 );
  WStack = bsxfun(@minus, WStack, t );
  
  if ~hasAnyNan
    % Tomasi Kanade without/with the metric constraint
    % Affine camera matrix, MLE estimation (Tomasi Kanade)
    % Reference: HZ2, p437, Algorithm 18.1
    
    [ U S V ] = svd( WStack );
    P = permute(reshape(U(:,1:3),2,nFrame,3),[1 3 2]);
    
    S=bsxfun(@times,[ S(1,1); S(2,2); S(3,3) ], V(:,1:3)');
  else
    [PStack,S]=lowRankDecomposition(WStack,3);
    
    P(1:2,1:3,:)=permute(reshape(PStack,2,nFrame,3),[1,3,2]);
  end
  % assign the translation
  P(3,4,:) = 1; P(1:2,4,:) = reshape( t, 2, 1, nFrame );
end
