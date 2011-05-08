function [ W, P, anim ]=generateW( anim, varargin )
% Generate projected data
%
% It uses the camera info and the 3D data to project to 2D
%
% USAGE
%  anim = anim.generateW()
%
% INPUTS
%  anim     - Animation object (help Animation for details)
%  varargin   - list of paramaters in quotes alternating with their values
%    P        - projection matrix. If not given, anim.P is used, if it does
%               not exist, it is computed from K,R,t
%    doFillP  - if set to true and anim.P is empty, it will fill anim.P
%    doFillW  - if set to true, will overwrite/create anim.
%
% OUTPUTS
%  W        - projected shape
%  P        - projection matrix
%  anim     - modified Animation
%
% EXAMPLE
%
% See also GENERATETOYANIMATION
%
% Vincent's Structure From Motion Toolbox      Version 3.1
% Copyright (C) 2008-2011 Vincent Rabaud.  [vrabaud-at-cs.ucsd.edu]
% Please email me if you find bugs, or have suggestions or questions!
% Licensed under the GPL [see external/gpl.txt]

[ P, doFillP, doFillW ] = getPrmDflt( varargin, ...
  { 'P' [] 'doFillP' false 'doFillW' false }, 1 );

W=[];

if isempty(anim.S)
  warning('S is not defined nor computable'); return;
end

if isempty(P)
  if isempty(anim.P); [ P anim ]=generateP(anim, doFillP);
  else P=anim.P;
  end
end

if ~isempty(anim.S)
  if size(P,2)==3
    % homography case
    W = bsxfun(@plus, multiTimes( P(:,1:2,:), anim.S, 1 ), P(:,3,:));
    W = reshape(normalizePoint(reshape(W,3,[]),3),2,anim.nPoint,...
      anim.nFrame);
  else
    % compute the projected points for normal camera and 3D points
    if size(anim.S,3)==1
      W = bsxfun( @plus, multiTimes( P(:,1:3,:), anim.S, 1 ), P(:,4,:) );
    else
      W = bsxfun( @plus, multiTimes( P(:,1:3,:), anim.S, 2 ), P(:,4,:) );
    end
    if anim.isProj; W = reshape(normalizePoint(reshape(W,3,[]),3),2,...
        anim.nPoint,anim.nFrame);
    else W = W(1:2,:,:);
    end
  end
end

% set missing values to NaN
if ~isempty(anim.mask); W(:,find(~anim.mask(:)))=NaN; end

if doFillW; anim.W=W; end
if doFillP; anim.P=P; end
