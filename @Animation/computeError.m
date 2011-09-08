function [ err errFrame errTot SAbsoluteT anim ] = computeError(anim, ...
  varargin)
% Compute several SFM errors
%
% It can compute:
%  - the reprojection error
%  - the 3D error between several recovered 3D shapes and a ground truth 3D
%  shape
%
% USAGE
%  [ err errFrame errTot ] = computeError( 'reproj', isProj, varargin )
%  [ err errFrame errTot SBest ] = computeError( '3D', isProj, varargin)
%
% INPUTS
%  anim       - Animation object
%  varargin   - list of paramaters in quotes alternating with their values
%       - 'doFillP', false, will replace/create anim.P with the projection
%       matrices
%       - 'animGT'   [] the ground truth anim to compare to
%       - 'doCheckAmbiguity' [true] if true, check for ambiguities due to
%       the camera model (scale for isProj and necker reversal for ~isProj)
%       - 'checkTransform' ['camera']. Transform that can be applied to get
%       a better match. Can be 'rigid', 'rigid+scale', 'homography' or
%       'camera' (cf the alignTo function)
%
% OUTPUTS 'reproj'
%  err       - [ 1 x 3 ] array of errors:
%              if no animGT is given:
%              err = mean( errFrame, 2 );
%              if an animGT is givne:
%              err(1) = mean( errFrame{1}(3,:), 2 )
%              err(2) = mean( errFrame{2}(3,:), 2 )
%              err(3) = mean( errFrame{3} )
%  errFrame  - if no animGT is given:
%              errFrame(1,:) = sum of squared reprojected distance
%                              differences per frame
%              errFrame(2,:) = sum of reprojected distances per frame
%                              divided by the span of the object
%              errFrame(3,:) = same as first one divided by the sum of
%                              squared distances
%  errFrame  - if an animGT is given, cell object:
%              {1} [ 3 x nFrame ] Bregler, CVPR00 error (SSD)
%              The rows contain the sums over the features (and for each
%              frame) of the squared errors in X-Y, in Z, and in 3D.
%              I.e. the rows contain for each frame:
%                sum_{i=1}^nPoint  (xReconstructed_i-xOriginal_i)^2
%                                         +(yReconstructed_i-yOriginal_i)^2
%                sum_{i=1}^nPoint  (zReconstructed_i-zOriginal_i)^2
%                row1+row2
%              {2} [ 3 x nFrame ] Torresani, Bregler, CVPR01 error
%              Same as above but distances are considered (not squared
%              distances), average and not sum, and the results are
%              normalized by the span of the object
%              {3} [ 1 x nFrame ] Torresani, PAMI08 error
%              At each frame:
%              \Vert SReconstructed-SOriginal \Vert_F/
%                               \Vert SOriginal \Vert_F
%  errTot   - [3 x nPoint x nFrame] error for each point at each frame
%             as a distance (simply abs of diff, not squared)
%  SAbsoluteT - the best transformed S to match SGTAbsolute (SGT
%               transformed by its R and t)
%  anim     - the original Animation object but modified is aksed to be
%
% EXAMPLE
%
% See also Animation
%
% Vincent's Structure From Motion Toolbox      Version 3.1.1
% Copyright (C) 2008-2011 Vincent Rabaud.  [vrabaud-at-cs.ucsd.edu]
% Please email me if you find bugs, or have suggestions or questions!
% Licensed under the GPL [see external/gpl.txt]

global IS_EXIST_PDIST;
if isempty(IS_EXIST_PDIST); IS_EXIST_PDIST=exist('pdist','file'); end

[ doFillP animGT doCheckAmbiguity checkTransform ] = ...
  getPrmDflt( varargin, ...
  { 'doFillP' false 'animGT' [] ...
  'doCheckAmbiguity' true 'checkTransform' 'camera'}, 1 );

nFrame = anim.nFrame; nPoint = anim.nPoint;
method='reproj';
if ~isempty(animGT); method = '3D'; end

switch method
  case 'reproj'
    % compute the projected points
    [ W disc anim ]=generateW(anim,'doFillP', doFillP);
    
    % precompute stuff for error number 2
    spanW = zeros(1,nFrame);
    if IS_EXIST_PDIST
      for i = 1 : nFrame
        spanW(i) = max(pdist( anim.W(:,:,i)));
      end
    else
      for i = 1 : nFrame
        spanW(i) = max(max(pdist2( anim.W(:,:,i), anim.W(:,:,i))));
      end
    end
    
    % compute the different errors
    errFrame = zeros( 3, nFrame );
    tmp = W - anim.W;
    
    % if there is a mask, only keep certain points on certain frames
    WNorm = W;
    if ~isempty(anim.mask)
      % set to 0 the ones that are 0 in the mask
      mask = repmat(reshape(~anim.mask,1,nPoint,nFrame), [2, 1, 1]);
      tmp(mask) = 0;
      WNorm(mask) = 0;
    end
    WNorm = sum(reshape(W,[],nFrame).^2, 1);
    
    errTot = abs(tmp);
    errFrame(1,:) = sum(reshape(tmp,[],nFrame).^2,1);
    errFrame(2,:) = sum(reshape(errTot,[],nFrame),1)./spanW;
    errFrame(3,:) = sqrt(errFrame(1,:)./WNorm);
    
    err = mean( errFrame, 2 );
  case '3D'
    [ disc disc disc disc SAbsoluteT SGTAbsolute ] = alignTo(anim, ...
      animGT, checkTransform);
    
    SDiff=SAbsoluteT-SGTAbsolute;
    
    % Update the errors
    errTot=SDiff;
    
    errFrame=cell(1,3);
    % compute the first type of error
    %            {1} [ 3 x nFrame ] Bregler, CVPR00 error (SSD)
    %            The rows contain the sums over the features (and for each
    %            frame) of the squared errors in X-Y, in Z, and in 3D.
    %            I.e. the rows contain for each frame:
    %              sum_{i=1}^nPoint  (xReconstructed_i-xOriginal_i)^2
    %                                     +(yReconstructed_i-yOriginal_i)^2
    %              sum_{i=1}^nPoint  (zReconstructed_i-zOriginal_i)^2
    %                row1+row2
    errFrame{1} = zeros( 3, nFrame );
    errFrame{1}(1,:) = reshape(sum(sum( SDiff(1:2,:,:).^2, 2 ),1),1,[]);
    errFrame{1}(2,:) = reshape(sum(SDiff(3,:,:).^2, 2 ),1,[]);
    errFrame{1}(3,:) = errFrame{1}(1,:)+errFrame{1}(2,:);
    
    %       {2} [ 3 x nFrame ] Torresani, Bregler, CVPR01 error
    %       Same as above but distances are considered (not squared
    %       distances), average and not sum, and the results are normalized
    %       by the span of the object
    normS = sqrt(reshape(sum(sum(animGT.S.^2,1),2),1,[]));
    errFrame{2} = zeros( 3, nFrame );
    errFrame{2}(1,:) = mean(sqrt(reshape(sum(sum( SDiff(1:2,:,:).^2, ...
      2),1),1,[])));
    errFrame{2}(2,:) = mean(reshape(SDiff(3,:,:),1,[]));
    errFrame{2}(3,:) = mean(sqrt(reshape(sum(sum( SDiff(:,:,:).^2, ...
      2),1),1,[])));
    % precompute stuff for error number 2
    spanS = zeros(1,size(animGT.S,3));
    if exist('pdist','file')
      for i = 1 : size(animGT.S,3)
        spanS(i) = max(pdist( animGT.S(:,:,i)));
      end
    else
      for i = 1 : size(animGT.S,3)
        spanS(i) = max(max(pdist2( animGT.S(:,:,i), animGT.S(:,:,i))));
      end
    end
    if size(spanS,2)==1; spanS = repmat(spanS,1,nFrame); end
    errFrame{2} = bsxfun(@rdivide, errFrame{2}, spanS);
    
    %              {3} [ 1 x nFrame ] Torresani, PAMI08 error
    %              At each frame:
    %              \Vert SReconstructed-SOriginal \Vert_F/
    %                               \Vert SOriginal \Vert_F
    errFrame{3} = sqrt(errFrame{1}(3,i))./normS;
    
    % Update other errors
    err=[ mean(errFrame{1}(3,:)) mean(errFrame{2}(3,:)) ...
      mean(errFrame{3}) ];
end
