function [ R, t, s, chirality, SAbsoluteT, SGTAbsolute ] = alignTo( ...
  anim, animGT, checkTransform)
% Align an animation to another by removing camera ambiguities
%
% For a projective camera, this function finds the best scale to multiply
% anim by to get it as close to animGT as possible.
% For an orthographic camera, this function finds the best depth and deals
% with the chirality issue, to get anim as close to animGT as possible.
%
% Basically, with the notations below:
%   SAbsoluteT = chirality*s*R*SAbsolute+t, SGTAbsolute = RGT*SGT + TT
% and SAbsoluteT and SGTAbsolute are as close to each other as possible
% chirality is a chirality matrix: diag([1 1 +-1])
% SAbsolute is S already transformed by its own rotation/translation
%
% USAGE
%  [ R, t, s, chirality, SAbsoluteT ] = anim.alignTo( animGT )
%
% INPUTS
%  animGT         - Animation to align to (ground truth)
%  checkTransform - the kind of transform that can be applied to anim
%                   ['none']. Transform that can be applied to get a
%                   better match. Can be 'homography' or 'rigid'
%                   'rigid+scale' or 'camera' (in which case it only looks
%                   for scale ambiguity with a projective camera, and for
%                   depths/chirality issues for an orthographic camera)
%
% OUTPUTS
%  R           - rotation to apply to the original S
%  t           - translation to apply as defined above
%  s           - scale to apply
%  chirality   - [ 1 x nFrame ] boolean, yes if z has to be opposed
%  SAbsoluteT  - transformed S
%  SGTAbsolute - as defined above
%
% EXAMPLE
%
% See also
%
% Vincent's Structure From Motion Toolbox      Version 3.1
% Copyright (C) 2008-2011 Vincent Rabaud.  [vrabaud-at-cs.ucsd.edu]
% Please email me if you find bugs, or have suggestions or questions!
% Licensed under the GPL [see external/gpl.txt]

% compute the globally transformed 3D shapes
nFrameList = [];
if anim.nFrame > 1; nFrameList = [ anim.nFrame ]; end
if animGT.nFrame > 1; nFrameList(end+1) = animGT.nFrame; end

nFrame = min(nFrameList);
if ~isempty(anim.t) && ~isempty(anim.R)
  SAbsolute = generateSAbsolute(anim);
else
  SAbsolute = repmat(anim.S,[1 1 nFrame]);
end
if ~isempty(animGT.t) && ~isempty(animGT.R)
  SGTAbsolute = generateSAbsolute(animGT);
else
  if size(animGT.S,3)==nFrame
    SGTAbsolute = animGT.S;
  else
    SGTAbsolute = repmat(animGT.S,[1 1 nFrame]);
  end
end

% prepare the data to send back
R=[]; t=[]; s=[]; chirality=[];
switch checkTransform
  case 'homography',
    nCoord = size(SAbsolute,1);
    R=zeros(nCoord+1,nCoord+1,anim.nFrame);
  case 'rigid'
    R=zeros(3,3,anim.nFrame); t=zeros(3,anim.nFrame);
  case 'rigid+scale'
    R=zeros(3,3,anim.nFrame); t=zeros(3,anim.nFrame);
    s=zeros(1,anim.nFrame);
  case 'camera'
    if anim.isProj
      % projective only has a scale ambiguity s*S=SGT
      % minimizing the square difference same as solving
      % Si*(s*Si-SGTi)=0 or s*sum(S(:).^2)-sum(S.*SGT)=0
      s=reshape(sum(sum(SAbsolute.*SGTAbsolute,1),2)./ ...
        sum(sum(SAbsolute.^2,1),2), 1, anim.nFrame);
    else
      tChir=cell(1,2);
      meanS = mean(SAbsolute,2);
      tChir{2}=reshape(-meanS+mean(SGTAbsolute,2),3, ...
        anim.nFrame);
      meanS(3,:)=-meanS(3,:);
      tChir{1}=reshape(-meanS+mean(SGTAbsolute,2),3, ...
        anim.nFrame);
      % still need to check for chirality for each frame
      SAbsoluteT=zeros(3,anim.nPoint,anim.nFrame);
      chirality=false(1,anim.nFrame);
    end
end

% compute the best alignment depending on the case
for i = 1 : nFrame
  S = SAbsolute(:,:,i); SGT = SGTAbsolute(:,:,i);
  
  switch checkTransform
    case 'homography',
      R(:,:,i) = computeHomography( S, SGT, 1 );
    case 'rigid'
      [ R(:,:,i) t(:,i) ]=computeOrientation(SGT,S,...
        'absolute');
    case 'rigid+scale'
      % include scale
      [ R(:,:,i) t(:,i) s(i) ]=computeOrientation(SGT,S,...
        'absolute');
    case 'camera'
      if ~anim.isProj
        % everything has been computed, we just need to check for chirality
        errBest = Inf;
        for j=1:2
          STmp=S;
          if j==1; STmp(3,:)=-STmp(3,:); end
          STmp=bsxfun(@plus,STmp,tChir{j}(:,i));
          
          errTmp=norm(STmp-SGT);
          if errTmp<errBest
            errBest=errTmp; SAbsoluteT(:,:,i)=STmp;
            if j==1; chirality(i)=true; else chirality(i)=false; end
          end
        end
      end
  end
end

if nargout <=4; return; end

% coompute the transformed S
switch checkTransform
  case 'homography',
    SAbsoluteT=multiTimes(R,normalizePoint(SAbsolute,-(nCoord+1)),2);
    SAbsoluteT=normalizePoint(SAbsoluteT,nCoord+1);
  case 'rigid'
    SAbsoluteT=multiTimes(R,SAbsolute,2);
    SAbsoluteT=bsxfun(@plus,SAbsoluteT,reshape(t,3,1,nFrame));
  case 'rigid+scale'
    SAbsoluteT=bsxfun(@times,reshape(s,1,1,nFrame),...
      multiTimes(R,SAbsolute,2));
    SAbsoluteT=bsxfun(@plus,SAbsoluteT,reshape(t,3,1,nFrame));
  case 'camera'
    if anim.isProj
      SAbsoluteT=bsxfun(@times,reshape(s,1,1,nFrame),SAbsolute);
    end
end
