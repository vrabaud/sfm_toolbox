function varargout=subsref(anim,idx)
% Returns a a specific member of the class
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
% Copyright (C) 2009 Vincent Rabaud.  [vrabaud-at-cs.ucsd.edu]
% Please email me if you find bugs, or have suggestions or questions!
% Licensed under the GPL [see external/gpl.txt]

if length(idx)==1
  % measurements
  switch idx.subs
    case 'W'
      varargout{1}=anim.W;
    case 'mask'
      varargout{1}=anim.mask;
      % 3D / Object
    case 'S'
      varargout{1}=anim.S;
    case 'conn'
      varargout{1}=anim.conn; % Cell of arrays of connectivities
      % NRSFM
    case 'l'
      varargout{1}=anim.l;
    case 'SBasis'
      varargout{1}=anim.SBasis;
      % Camera
    case 'P'
      varargout{1}=anim.P;
    case 'K'
      varargout{1}=anim.K;
    case 'R'
      varargout{1}=anim.R;
    case 't'
      varargout{1}=anim.t;
      % Misc
    case 'misc'
      varargout{1}=anim.misc;
    case 'isProj'
      % Info
      varargout{1}=anim.isProj;
      % Quantities
    case 'nBasis'
      varargout{1}=anim.nBasis;
    case 'nPoint'
      varargout{1}=anim.nPoint;
    case 'nFrame'
      varargout{1}=anim.nFrame;
  end
else
  % measurements
  switch idx(1).subs
    case 'addNoise',
      varargout{1}=addNoise(anim,idx(2).subs{:});
    case 'alignTo',
      switch nargout
        case 1,
          varargout{1} = alignTo(anim,idx(2).subs{:});
        case 2,
          [ varargout{1}, varargout{2} ] = alignTo(anim,idx(2).subs{:});
        case 3,
          [ varargout{1}, varargout{2}, varargout{3} ] = alignTo(...
            anim,idx(2).subs{:});
        case 4,
          [ varargout{1}, varargout{2}, varargout{3}, varargout{4} ] = ...
            alignTo(anim,idx(2).subs{:});
        case 5,
          [ varargout{1}, varargout{2}, varargout{3}, varargout{4}, ...
            varargout{5} ] = alignTo(anim,idx(2).subs{:});
        case 6,
          [ varargout{1}, varargout{2}, varargout{3}, varargout{4}, ...
            varargout{5}, varargout{6} ] = alignTo(anim,idx(2).subs{:});
      end
    case 'centerW',
      varargout{1}=centerW(anim);
    case 'computeError',
      switch nargout
        case 1,
          varargout{1} = computeError(anim,idx(2).subs{:});
        case 2,
          [ varargout{1}, varargout{2} ] = computeError(anim, ...
            idx(2).subs{:});
        case 3,
          [ varargout{1}, varargout{2}, varargout{3} ] = computeError(...
            anim,idx(2).subs{:});
        case 4,
          [ varargout{1}, varargout{2}, varargout{3}, varargout{4} ] = ...
            computeError(anim,idx(2).subs{:});
        case 5,
          [ varargout{1}, varargout{2}, varargout{3}, varargout{4}, ...
            varargout{5} ] = computeError(anim,idx(2).subs{:});
      end
    case 'generateCamFromRt',
      varargout{1}=generateCamFromRt(anim);
    case 'generateP',
      switch nargout
        case 1,
          varargout{1} = generateP(anim,idx(2).subs{:});
        case 2,
          [ varargout{1}, varargout{2} ] = generateP(anim,idx(2).subs{:});
      end
    case 'generateSAbsolute',
      varargout{1} = generateSAbsolute(anim);
    case 'generateW',
      switch nargout
        case 1,
          varargout{1} = generateW(anim,idx(2).subs{:});
        case 2,
          [ varargout{1}, varargout{2} ] = generateW(anim,idx(2).subs{:});
        case 3,
          [ varargout{1}, varargout{2}, varargout{3} ] = generateW(anim,...
            idx(2).subs{:});
      end
    case 'sampleFrame',
      switch nargout
        case 1,
          varargout{1} = sampleFrame(anim,idx(2).subs{:});
        case 2,
          [ varargout{1}, varargout{2} ] = sampleFrame(anim, ...
            idx(2).subs{:});
      end
    case 'setFirstRToId',
      varargout{1}=setFirstRToId(anim);
    case 'W'
      varargout{1}=subsref(anim.W,idx(2));
    case 'mask'
      varargout{1}=subsref(anim.mask,idx(2));
      % 3D / Object
    case 'S'
      varargout{1}=subsref(anim.S,idx(2));
    case 'conn'
      % Cell of arrays of connectivities
      varargout{1}=subsref(anim.conn,idx(2));
      % NRSFM
    case 'l'
      varargout{1}=subsref(anim.l,idx(2));
    case 'SBasis'
      varargout{1}=subsref(anim.SBasis,idx(2));
      % Camera
    case 'P'
      varargout{1}=subsref(anim.P,idx(2));
    case 'K'
      varargout{1}=subsref(anim.K,idx(2));
    case 'R'
      varargout{1}=subsref(anim.R,idx(2));
    case 't'
      varargout{1}=subsref(anim.t,idx(2));
      % Misc
    case 'misc'
      varargout{1}=subsref(anim.misc,idx(2));
  end
end
