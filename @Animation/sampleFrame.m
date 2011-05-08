function [ anim SimMin ] = sampleFrame( anim, nSamp, SimMin )
% Sample an Animation by reducing the number of frames
%
% Given an Animation, it returns an Animation with nSamp frames
% The first and the last frame are kept, the other ones are as far
% appart from each other.
% If an array is given, it is assume to be the range to choose from
%
% If a SimMin matrix is given, it is also sampled inthe same way
%
% USAGE
%  anim = anim.sampleFrame( nSamp, SimMin )
%
% INPUTS
%  anim     - Animation object (help Animation for details)
%  nSamp    - number of frames to keep after sampling or range
%  SimMin   - [ anim.nFrame x anim.nFrame ] matrix that needs to
%             be sampled too
%
% OUTPUTS
%  anim     - modified Animation object
%  SimMin   - [ nSamp x nSamp ] sampled SimMin
%
% EXAMPLE
%
% See also GENERATETOYANIMATION
%
% Vincent's Structure From Motion Toolbox      Version 3.0
% Copyright (C) 2008-2011 Vincent Rabaud.  [vrabaud-at-cs.ucsd.edu]
% Please email me if you find bugs, or have suggestions or questions!
% Licensed under the GPL [see external/gpl.txt]

if length(nSamp)==1
  samp=round(1:(anim.nFrame-1)/(nSamp-1):anim.nFrame);
else
  samp = nSamp;
end

% figure out if P has to be resized
doP=isempty(anim.t) || isempty(anim.R);

s=fieldnames(anim);
for i=1:length(s)
  if ~isempty(anim.(s{i})) && ~isempty(anim.(s{i})) && ...
      ~strcmp(s{i},'SBasis')
    if ndims(anim.(s{i}))==3 && (~strcmp(s{i},'P') || doP)
      anim=subsasgn(anim,struct('type','.','subs',s{i}),...
        anim.(s{i})(:,:,samp));
    end
    if ndims(anim.(s{i}))==2 && ~strcmp(s{i},'S') && size(anim.(s{i}),2)>1
      anim=subsasgn(anim,struct('type','.','subs',s{i}),...
        anim.(s{i})(:,samp));
    end
  end
end

if nargin>=3
  SimMin = SimMin(samp,samp);
end
