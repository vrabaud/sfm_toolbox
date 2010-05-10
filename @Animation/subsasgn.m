function anim=subsasgn(anim,idx,rhs)
% Assigns a specific member of the class
%
% USAGE
%  anim = Animation()
%
% INPUTS
%
% OUTPUTS
%  anim     - an Animation object
%
% EXAMPLE
%
% Vincent's Structure From Motion Toolbox      Version 3.0
% Copyright (C) 2008-2010 Vincent Rabaud.  [vrabaud-at-cs.ucsd.edu]
% Please email me if you find bugs, or have suggestions or questions!
% Licensed under the GPL [see external/gpl.txt]

if max(size(idx))==1
  var=idx.subs;
  % measurements
  switch var
    case 'W'
      anim.W=rhs;
    case 'mask'
      anim.mask=logical(mask);
      % 3D / Object
    case 'S'
      % do not modify if shape basis is defined
      if ~isempty(anim.l) && ~isempty(anim.SBasis)
        warning('Cannot modify S as l and SBasis define it'); return;
      end
      anim.S=rhs;
    case 'conn'
      anim.conn=rhs; % Cell of arrays of connectivities
      % NRSFM
    case 'l'
      anim.l=rhs;
    case 'SBasis'
      anim.SBasis=rhs;
      % Camera
    case 'P'
      % do not modify if R and t are defined
      if ~isempty(anim.R) && ~isempty(anim.t)
        warning('Cannot modify P as R and t define it'); return;
      end
      anim.P=rhs;
    case 'K'
      anim.K=rhs;
    case 'R'
      anim.R=rhs;
    case 't'
      anim.t=rhs;
      % Misc
    case 'misc'
      anim.misc=rhs;
      % Info
    case 'isProj'
      anim.isProj=rhs;
  end
else
  var=idx(1).subs;
  % measurements
  switch var
    case 'W'
      anim.W=subsasgn(anim.W,idx(2),rhs);
    case 'mask'
      anim.mask=logical(subsasgn(anim.mask,idx(2),rhs));
      % 3D / Object
    case 'S'
      if ~isempty(anim.l) && ~isempty(anim.SBasis)
        warning('Cannot modify S as l and SBasis define it'); return;
      end
      anim.S=subsasgn(anim.S,idx(2),rhs);
    case 'conn'
      % Cell of arrays of connectivities
      anim.conn=subsasgn(anim.conn,idx(2),rhs);
      % NRSFM
    case 'l'
      anim.l=subsasgn(anim.l,idx(2),rhs);
    case 'SBasis'
      anim.SBasis=subsasgn(anim.SBasis,idx(2),rhs);
      % Camera
    case 'P'
      % do not modify if R and t are defined
      if ~isempty(anim.R) && ~isempty(anim.t)
        warning('Cannot modify P as R and t define it'); return;
      end
      anim.P=subsasgn(anim.P,idx(2),rhs);
    case 'K'
      anim.K=subsasgn(anim.K,idx(2),rhs);
    case 'R'
      anim.R=subsasgn(anim.R,idx(2),rhs);
    case 't'
      anim.t=subsasgn(anim.t,idx(2),rhs);
      % Misc
    case 'misc'
      anim.misc=subsasgn(anim.misc,idx(2),rhs);
  end
end

% modify (or set to be modified) interdependent elements
switch var
  % deal with the modified shape basis/coefficients
  case {'l','SBasis'}
    anim.nFrame=max([ size(anim.W,3) size(anim.S,3) size(anim.l,2) ...
      size(anim.R,3) ] );
    if ~isempty(anim.l) && ~isempty(anim.SBasis)
      anim.nBasis=size(anim.SBasis,3);
      anim=generateSFromLSBasis(anim);
    end
  case {'K', 'R', 't'}
    % reset some global values
    if size(anim.K,1)==5; anim.isProj=true;
    elseif size(anim.K,1)==3; anim.isProj=false;
    end
    anim.nFrame=max([ size(anim.W,3) size(anim.S,3) size(anim.l,2) ...
      size(anim.R,3) ] );
    % yeah yeah, not optimal but easier to keep everything in sync like
    % that
    if ~isempty(anim.R) && size(anim.t,2)==size(anim.R,3)
      anim.P=generateP(anim);
    end
  case {'W', 'S'}
    anim.nPoint=max([ size(anim.W,2) size(anim.S,2) size(anim.SBasis,2)] );
    anim.nFrame=max([ size(anim.W,3) size(anim.S,3) size(anim.l,2) ...
      size(anim.R,3) ] );
end
