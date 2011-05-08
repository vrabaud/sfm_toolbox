function display(anim)
% Display the different elements contained in the object
%
% USAGE
%  anim.display()
%
% INPUTS
%
% OUTPUTS
%  anim     - an Animation object
%
% EXAMPLE
%
% Vincent's Structure From Motion Toolbox      Version 3.0
% Copyright (C) 2008-2011 Vincent Rabaud.  [vrabaud-at-cs.ucsd.edu]
% Please email me if you find bugs, or have suggestions or questions!
% Licensed under the GPL [see external/gpl.txt]

fprintf('%s is an Animation object from Vincent''s SFM toolbox\n', ...
  inputname(1))

% measurements
displaySize('W',size(anim.W))
displaySize('mask',size(anim.mask))

% 3D / Object
displaySize('S',size(anim.S))
displaySize('conn',size(anim.conn))
% NRSFM
displaySize('l',size(anim.l))
displaySize('SBasis',size(anim.SBasis))
% Camera
displaySize('P',size(anim.P))
displaySize('K',size(anim.K))
displaySize('R',size(anim.R))
displaySize('t',size(anim.t))
% Misc
displaySize('misc',size(anim.misc))
% Info
displayInt('isProj',anim.isProj)
displayInt('type',anim.type)
% Quantities
displayInt('nBasis',anim.nBasis)
displayInt('nPoint',anim.nPoint)
displayInt('nFrame',anim.nFrame)

  function displaySize(s,siz)
    switch length(siz)
      case 2,
        if min(siz)==0
          fprintf('%s : []\n',s)
        else
          fprintf('%s : %i x %i\n',s,siz(1),siz(2))
        end
      case 3,
        fprintf('%s : %i x %i x %i\n',s,siz(1),siz(2),siz(3))
    end
  end
  function displayInt(s,val)
    fprintf('%s : %i\n',s,val)
  end
end
