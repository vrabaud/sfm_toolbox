function [arg1, arg2]=extractFromTFT(T,type)
% Extract information from the trifocal tensor
%
% USAGE
%  [ep, epp]=extractFromTFT(T)
%
% INPUTS 1 - returns the epipoles
%  T       - 3x3x3 trifocal tensor
%  type    - 0
%
% OUTPUTS 1
%  ep      - epipole second picture
%  epp     - epipole third picture
%
% EXAMPLE
%
% See also
%
% Vincent's Structure From Motion Toolbox      Version 1.1
% Copyright (C) 2008-2011 Vincent Rabaud.  [vrabaud-at-cs.ucsd.edu]
% Please email me if you find bugs, or have suggestions or questions!
% Licensed under the GPL [see external/gpl.txt]

switch type
  case 0
    % Extract the epipoles of the 2 other views
    % Reference: HZ2, p.395
    V=zeros(3,3); for i=1:3; V(i,:)=solveLeastSqAx(T(:,:,i))'; end
    epp=solveLeastSqAx(V);
    for i=1:3; V(i,:)=solveLeastSqAx(T(:,:,i)')'; end
    ep=solveLeastSqAx(V);
    
    arg1=ep; arg2=epp;
end
